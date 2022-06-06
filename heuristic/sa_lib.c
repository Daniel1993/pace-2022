#include <math.h>
#include <assert.h>

#include "sa_lib.h"
#include "aux.h"
#include "WELL512a.h"
#include "simulatedAnnealing.h"

SA_s internal_lib_sa_;
FILE* internal_prof_file_;
double thresholdOpsNotImproving = THRESHOLD_OPS_NOT_IMPROVING; // increase as time passes
double thresholdResetTemperature = THRESHOLD_RESET_TEMPERATURE;

/*extern */float reorganize_factor = MAX_REORGANIZE_FACTOR;

static const double NB_BATCHES_IN_EXEC =
#if defined(BATCH_SIZE) && BATCH_SIZE > 0
  (double)BATCH_SIZE
#else
  (double)1.0
#endif
;

int
SA_get_nbVerts()
{
  SA_s sa = internal_lib_sa_;
  return sa->G->v;
}

int
SA_get_nbEdges()
{
  SA_s sa = internal_lib_sa_;
  return sa->G->e;
}

int
SA_init_G(SA_parameters_s params, graph G)
{
  SA_s sa = (SA_s)calloc(1, sizeof(struct SA_));
  internal_lib_sa_ = sa;
  internal_prof_file_ = stdout;

  // RNG initialization
  uint32_t init_state[512/(sizeof(uint32_t)*8)]; // 512 bits
  WriteEntropyInBuffer(init_state, 512/8);
  InitWELLRNG512a(init_state);

  sa->hot = params.hot;
  sa->hotD = params.hotD;
  sa->cold = params.cold;
  sa->coldD = params.coldD;
  sa->decay = params.decay;

  G->t = NULL;
  G->d = NULL;
  sa->G = G;
  sa->order = malloc((1+sa->G->v)*sizeof(int));
  sa->orderTransSCC = malloc(sa->G->v*sizeof(int));
  sa->SCCsize = malloc((1+sa->G->v)*sizeof(int));
  sa->SCCnbEdges = malloc((1+sa->G->v)*sizeof(int));
  
  // fprintf(stderr, "Before sccRestrict graph (v=%li,e=%li)\n", G->v, G->e);
  sa->rG = sccRestrict(sa->G, sa->order, sa->SCCsize, sa->SCCnbEdges, 0);
  // fprintf(stderr, "Before contraction graph (#SCCs=%i)\n", sa->order[0]);
  // TODO: use Levy and Low to contract the graph
  int rep = 1, nbEdges, oldFVSs;
  while (rep)
  {
    oldFVSs = sa->rG->inFVSs;
    nbEdges = sa->rG->e;
    rep = contractG(sa->rG);
    // fprintf(stderr, "After contractG: new FVSsize = %d (old edges=%d, new edges=%d)\n", sa->rG->inFVSs, nbEdges, sa->rG->e);
    if (sa->rG->inFVSs > oldFVSs)
    { // only rerun sccRestrict in case there are removes that can cause SCC splits
      nbEdges = sa->rG->e;
      sccRestrict(G, sa->order, sa->SCCsize, sa->SCCnbEdges, 1);
      if (nbEdges != sa->rG->e)
        rep = 1; // could be there is splits
    }
  }

  sa->SCCgraphs = (graph*)malloc(sa->order[0]*sizeof(graph));

  sa->maxE = calloc(sa->order[0], sizeof(int));
  sa->localMaxE = calloc(sa->order[0], sizeof(int));
  sa->currMaxE = calloc(sa->order[0], sizeof(int));
  sa->s = calloc(sa->order[0], sizeof(sa_state));
  sa->vertsPerSCC  = calloc(sa->order[0]+1, sizeof(int*));
  // the other SCC pointers point to the same array but in different locations
  sa->vertsPerSCC[0] = malloc((sa->G->v*2)*sizeof(int));

  sa->complexSCCs = (int*)malloc((sa->order[0]+1)*sizeof(int));
  int *complexSCCs = sa->complexSCCs;

  int *A = &(sa->order[1]);
  int *At = &(sa->orderTransSCC[0]);
  sa->rG->d = (int*)malloc(sa->rG->v*sizeof(int)); // TODO: move this elsewhere
  for (int i = 0; i < sa->order[0]; ++i)
  {
    // fprintf(stderr, "Before extractSCC%i\n", i);
    sa->currSCC = i;
    sa->SCCgraphs[i] = extractSCC(sa->rG, A, sa->SCCsize[i]);
    // fprintf(stderr, "After extractSCC%i (v=%i,e=%i)\n", i, sa->SCCgraphs[i]->v, sa->SCCgraphs[i]->e);
    sa->s[i] = allocS(sa->SCCgraphs[i]);
    for (int j = 0; j < sa->SCCsize[i]; ++j)
      At[j] = j; // ids start from 0 in those SCCs, use G->d[>scc_vert_id>] to translate
    prepareS_greedy(sa->s[i], At, sa->SCCsize[i]);
    A += sa->SCCsize[i];
    At += sa->SCCsize[i];
    if (i > 0)
      sa->vertsPerSCC[i] = sa->vertsPerSCC[i-1] + sa->SCCsize[i-1] + 1;
    
    if(sa->SCCnbEdges[i] == sa->SCCsize[i]*sa->SCCsize[i]-sa->SCCsize[i])
    {
      //! Deal with complete graphs
      *(sa->vertsPerSCC[i]) = 0; // local SCC graph vert index
      *(sa->vertsPerSCC[i]+1) = -1; // end of DAG
      sa->maxE[i] = 1;
      sa->localMaxE[i] = 1;
      sa->currMaxE[i] = 1;
    }
    else
    {
      *(complexSCCs++) = i;
      guessSA(sa->s[i], sa->vertsPerSCC[i], /* TODO: another paramenter */1);
      // shuffleSA(sa->s[i], sa->SCCnbEdges, sa->SCCsize, i, sa->maxE[i]);
    }

  }
  *complexSCCs = -1;
  return 0;
}

int
SA_init_F(SA_parameters_s params, const char *filename)
{

  FILE *stream = NULL;
  if (filename)
    stream = fopen(filename, "r");
  if (!stream)
    stream = stdin;
  graph G = loadG(stream);
  if (stream != stdin)
    fclose(stream);

  return SA_init_G(params, G);
}

void
SA_destroy()
{
  SA_s sa = internal_lib_sa_;
  for (int i = 0; i < sa->order[0]; ++i) {
    freeS(sa->s[i]);
    freeG(sa->SCCgraphs[i]);
  }
  free(sa->maxE);
  free(sa->localMaxE);
  free(sa->currMaxE);
  free(sa->vertsPerSCC[0]);
  free(sa->vertsPerSCC);
  free(sa->s);
  freeG(sa->G);
  free(sa->order);
  free(sa->SCCsize);
  free(sa);
}

int*  // loop with while(-1 != *bestSolution)
SA_getBestSolution()
{
  SA_s sa = internal_lib_sa_;
  return getBestSolutionBuffer(sa->G);
}

void
SA_run()
{
  SA_s sa = internal_lib_sa_;
  double hotStart = getTemperature(sa->hot/(float)(sa->G->v), sa->hotD);
  double coldStart = getTemperature(sa->cold/(float)(sa->G->v), sa->coldD);

  double hot;

  int maxE = 0;

  int *L = getBestSolutionBuffer(sa->G); /* Current best config */
  int *P;

  sa->opsNotImproving = 0L;

  hot = coldStart; // the first run do with low temperature
  reorganize_factor = 1.2f;
  int firstTime = 1;
  while (1)
  {
    int *complexSCCs = sa->complexSCCs;
    long threshold = 50 + thresholdResetTemperature * sa->G->v;
    if (-1 == *complexSCCs)
      break; // no complex SCCs to process
    while (-1 != *complexSCCs)
    {
      int i = *(complexSCCs++);
      int *vertsPerSCC = sa->vertsPerSCC[i];

      sa->currSCC = i;
      sa->currSCCsize = sa->SCCsize[i];

#ifndef NDEBUG
      vertsPerSCC = 
#endif
      executeSA(
        sa->s[i],
        hot,
        sa->decay,
        vertsPerSCC,
        firstTime
      );
      assert(-1 == *vertsPerSCC && "Wrong list ending");
      assert(sa->currSCCsize > vertsPerSCC - sa->vertsPerSCC[i] && "Too many vertexes in list");
    } // end loop per SCC
    // // ----------------
    // int bestSol = 0;
    // int maxE = 0, maxLocE = 0, currE = 0;
    // for (int i = 0; i < sa->order[0]; i++) {
    //   int *vertsPerSCC = sa->vertsPerSCC[i];
    //   maxE += sa->maxE[i];
    //   maxLocE += sa->localMaxE[i];
    //   currE += sa->currMaxE[i];
    //   while (-1 != *vertsPerSCC)
    //   {
    //     vertsPerSCC++;
    //     bestSol++;
    //   }
    // }
    // fprintf(stderr, "temp=%e, ops=%li, factor=%f bestSol=%li E=%i lE=%i cE=%i\n",
    //   hot, sa->opsNotImproving, reorganize_factor, sa->G->v - bestSol + sa->rG->inFVSs, maxE, maxLocE, currE);
    // // ----------------
    if (sa->opsNotImproving < threshold)
    {
      if (!firstTime)
      {
        if (hot > coldStart)
          hot *= sa->decay;
        if (reorganize_factor > MIN_REORGANIZE_FACTOR)
          reorganize_factor *= REORGANIZE_FACTOR_DECAY;
      }
    }
    else
    {
      firstTime = 0;
      // fprintf(stderr, " --------- RESET_T --------- \n");
      int maxE = 0, currE = 0, maxLocE = 0;
      for (int i = 0; i < sa->order[0]; i++) {
        maxE += sa->maxE[i];
        currE += getE(sa->s[i]);
        maxLocE += sa->localMaxE[i];
      }
      complexSCCs = sa->complexSCCs;
      while (-1 != *complexSCCs)
      {
        int i = *(complexSCCs++);
        sa->localMaxE[i] = 0;
        shuffleSA(sa->s[i], sa->SCCnbEdges, sa->SCCsize, i, sa->maxE[i]);
      }
      reorganize_factor = MAX_REORGANIZE_FACTOR;
      hot = hotStart;
      
      sa->opsNotImproving = 0; // reset;
      if (thresholdResetTemperature < MAX_THRESHOLD_RESET_TEMPERATURE)
      {
        thresholdOpsNotImproving += THRESHOLD_OPS_NOT_IMPROVING; // TODO
        thresholdResetTemperature += THRESHOLD_RESET_TEMPERATURE;
      }
    }
    if (sigterm_called)
      break; // timeout
  } // end loop 1

EXITPOS_SA:

#ifndef NDEBUG
  for (int i = 0; i < sa->order[0]; i++) {
    int *vertsPerSCC = sa->vertsPerSCC[i];
    assert(-1 != *vertsPerSCC); // at least 1 per SCC
    while (*(++vertsPerSCC) != -1);
    assert(sa->SCCsize[i] >= vertsPerSCC - sa->vertsPerSCC[i]);
  }
#endif

  P = L;
  for (int i = 0; i < sa->order[0]; i++)
  {
    int *vertsPerSCC = sa->vertsPerSCC[i];
    while (*vertsPerSCC != -1)
    {
      assert(sa->G->v > L - P);
      *L = sa->SCCgraphs[i]->d[*vertsPerSCC];
      assert(-1 != *L);
#ifndef NDEBUG
      for (int *kL = P; kL != P; ++kL) assert(*kL != *L);
#endif
      L++;
      vertsPerSCC++;
    }
  }
  *L = -1;
  sa->maxDagSize = L - P;
  L = P;
}

void
SA_set_prof_file(FILE *fp)
{
  internal_prof_file_ = fp;
}

void SA_printFVS()
{
  SA_s sa = internal_lib_sa_;
  // int sol = 0;
  unsigned char inMaxDag[sa->G->v];
  memset(inMaxDag, 0, sa->G->v*sizeof(unsigned char));
  for(int *L = SA_getBestSolution(); -1 != *L; L++) {
    assert(*L >= 0 && *L < sa->G->v && "wrong solution");
    inMaxDag[*L] = 1;
    // sol++;
    // if (!(*L >= 0 && *L < sa->G->v)) fprintf(stderr, "wrong sol : %d\n", *L);
  }
  for (int i = 0; i < sa->rG->inFVSs; ++i)
    inMaxDag[sa->rG->inFVS[i]] = 0; // remove the known FVS
  // fprintf(stderr, "sol : %d\n", sol);
  for(int i = 0; sa->G->v > i; i++)
    if (!inMaxDag[i])
      fprintf(internal_prof_file_, "%d\n", i+1);
}
