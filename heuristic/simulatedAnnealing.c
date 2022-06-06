#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fenv.h>
#include <assert.h>
#include <limits.h>
#include <time.h>

// #include <bsd/stdlib.h>
#include "WELL512a.h"
#include "sa_state.h"
#include "sa_lib.h"
#include "aux.h"

// extern unsigned long nbProcVert;
// extern unsigned long nbDfsCalls;

extern SA_s internal_lib_sa_;

/* Used for defining calibration values */
double
getTemperature(double p, /* probability */
	       int dE /* Energy delta */
	       )
{
  assert(0.5 > p && "Invalid probability to define temp.");
  double res = dE / log2((1/p)-1);
  // printf("temp: %f\n", res);
  return res;
}

// /* Used for defining calibration values */
// double // TODO: dE == 1
// getProb(double T, int d)
// { // TODO: check the formula
//   double res = 1 /(2^{d/T} - 1);
//   // printf("temp: %f\n", res);
//   return res;
// }

extern int flag_did_reorganize;
extern int flag_reorganize_improve;

static int
improveSolution(sa_state s,
  int *L,
  int **R,
  int doGetSketch
      )
{
  SA_s sa = internal_lib_sa_;
  int sccIdx = sa->currSCC;
  int localMaxImprove = 0;
  sa->currMaxE[sccIdx] = getE(s);

  // fprintf(stderr, "i");
  // do not take sketches on first time in SA
  if ( sa->currMaxE[sccIdx] > sa->localMaxE[sccIdx] )
  {
    localMaxImprove = 1;
    sa->localMaxE[sccIdx] = sa->currMaxE[sccIdx];
    if (sa->currMaxE[sccIdx] > sa->maxE[sccIdx])
    {
      sa->maxE[sccIdx] = sa->currMaxE[sccIdx];
      /* Store if you are about to lose best config. */
      if (doGetSketch)
      {
        *R = getSketch(s, L); // TODO: if (0 > pdE) ?
        // prepareS_empty(s);
      }
    }
  }
  return localMaxImprove;
}

static int
handleTemperature(double T
      )
{
  int dE  = INT_MAX; // Energy delta
  uint32_t B = WELLRNG512a_u4(); // replace for bsd function
  if (0 != B) 
  {
    float t; // temporary variable
    t = log2f(-B);
    t -= log2f(B);
    t *= T;
    dE = lrintf(t);
  }
  if(dE > 1)
    dE = 1; // Impossible to improve more than 1 in this problem
  return dE;
}

/* Return pointer to the end of the buffer */
int *
executeSA(
    sa_state s, /* State struct to use */
	  double hot, /* Hot temperature */
	  double decay, /* NOT USED */
	  int *L, /* Buffer for storing solution */
    int improveOnlyAfter /* if set to 1, only updates the final solution on exit */
	  )
{
  SA_s sa = internal_lib_sa_;
  int *R = L; /* Return value */
  double T = hot; /* Temperature */
  int pdE;
  int sccIdx = sa->currSCC;
  long threshold = 5 + thresholdOpsNotImproving * sa->SCCsize[sccIdx];
  long incReorganize = 1 + thresholdOpsNotImproving * sa->SCCsize[sccIdx];
  fesetround(FE_TONEAREST); /* For dE computation */
  int nbCandidates = getNbCandidates(s);
  int countImproves = 0;

  // function re-enters multiple times (compute temperature)
  for (long i = 0; ; i++, sa->opsNotImproving++)
  {
    int nbCandidates = getNbCandidates(s);
    for (long j = 0; j < nbCandidates; ++j)
    {
      int dE = handleTemperature(T);
      choose(s); /* Generate next sa_state */
      pdE = check(s, dE); /* Returns proposed dE */
      // if pdE == 1 accept else always reject --> Greedy
      if ( dE <= pdE )
      {
        if (improveSolution(s, L, &R, 1/* !improveOnlyAfter */))
          countImproves++;
        execute(s);
      }
      if ( sigterm_called )
        break;
    }
    choose(s); // restore candidate list in case it is empty

    if (countImproves)
      sa->opsNotImproving = 0;

    // fprintf(stderr, "Nb improvements: %d, ratio: %f\n", countImproves, (float)nbCandidates/(float)countImproves);
    countImproves = 0;

    if (!improveOnlyAfter && flag_did_reorganize && !flag_reorganize_improve)
    {
      i += incReorganize;
      sa->opsNotImproving += incReorganize;
    }

    tryReorganize(s);
    flag_did_reorganize = 0;
    flag_reorganize_improve = 0;
    if ( sigterm_called || i > threshold )
      goto EXITPOS_SA;

    /* if (improveOnlyAfter)
      shuffleCandidates(s); */
  }

EXITPOS_SA:
  if (improveOnlyAfter)
  { // TODO: this may make the solution decay
    R = getSketch(s, L);
    sa->maxE[sccIdx] = getE(s);
    sa->currMaxE[sccIdx] = sa->maxE[sccIdx];
    sa->localMaxE[sccIdx] = sa->maxE[sccIdx];
  }
#ifndef NDEBUG
  while (*R != -1)
    R++;
  assert(sa->currSCCsize > R - L && "Wrong sized solution");
#endif

  return R;
}

int *
guessSA(sa_state s,
  int *L,
  int attempts
      )
{
  SA_s sa = internal_lib_sa_;
  int *R = L; /* Return value */
  int sccIdx = sa->currSCC;

  while (0 < --attempts)
    decayAndReorganize(s);
  return R;
}

void
shuffleSA(
    sa_state s, /* State struct to use */
    int *SCCnbEdges,
    int *SCCsize,
    int curSCC,
    int nbVerts
	  )
{
  if (SCCnbEdges[curSCC] == SCCsize[curSCC]*SCCsize[curSCC]-SCCsize[curSCC])
    return;
  decayAndReorganize(s);
}
