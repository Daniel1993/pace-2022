#define _GNU_SOURCE // qsort_r()
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <math.h>

// #include <bsd/stdlib.h>
#include "WELL512a.h"
#include "graph.h"
#include "splayTree.h"
#include "darray.h"
#include "reorganize.h"
#include "sa_state.h"
#include "bitmap.h"

#ifdef POSTSCRIPT
FILE *LOG = NULL;
int POST_NB_VERT = 0;
#endif /* POSTSCRIPT */

#ifndef NDEBUG
static void
assertState(sa_state s);
#endif /* NDEBUG */

#define ARC4RAND_UNIF(_num) (WELLRNG512a_u4()%_num)
// #define ARC4RAND_UNIF(_num) arc4random_uniform(_num)
// #define ARC4RAND_UNIF(_num) random() % _num

// unsigned long nbProcVert = 0;

struct sa_state_ {
  graph G; /* The underlying graph */
  sTree t; /* K = vertex index. V = Topological Order index. */
  sTree O[2]; /* K = Topological index. V = sub counter. */
  int *Top2V; /* Maps topological order to vertex. */
  int *condemned; /* List of nodes to remove. Normaly small */
  int *scores; /* inDeg*outDeg */

  /* Candidate List stuff */
  // char *cBool; /* Candidate booleans */
  int cLs; /* Size of candidate list */
  int *cList; /* List of candidates */
  /* The current candidate is cList[0] */

  int *posInCLs; /* keeps position in candidate list */ 

  int cLs2; /* Check candidates in epochs */
  int *cList2;

  /* Excluded adjacency in check */
  darray DegEx[2]; /* Excluded vertexes */

  /* Position lives in the topological order */
  int p; /* The position selected by check */

  int c; // some counter for the check

  enum adjacencies dir; /* Direction selected by check */
  /* Delta E for last check. Used to trim execute. */
  int checkDE;
  /* int *inCList; */

#ifndef NDEBUG
  bitmap lastVerts;
  int sizeLastBestSketch;
#endif

  /* Supporting the reorganizing heuristic */
  organizer o;
};

sa_state
allocS(graph G)
{
  sa_state s = (sa_state) calloc(1, sizeof(struct sa_state_));

  s->G = G;
  s->t = allocSTree(1+G->v); /* Add a sentinel */
  s->O[in] = allocSTree(G->v); // TODO: these are used in the check state
  s->O[out] = allocSTree(G->v);
  s->Top2V = (int *)calloc(G->v, sizeof(int));
  s->condemned = (int *)calloc(2*G->v+1, sizeof(int));
  s->cList = (int *)malloc((G->v+1)*sizeof(int));
  s->cList2 = (int *)malloc((G->v+1)*sizeof(int));
  s->posInCLs = (int *)malloc(G->v*sizeof(int));
  memset(s->posInCLs, -1, G->v*sizeof(int));
  s->scores = (int *)malloc((2*G->v+1)*sizeof(int));
  for (int i = 0; i < G->v; ++i)
    if (G->deg[in][i] == 0 || G->deg[out][i] == 0)
      s->scores[i] = 0;
    else
      s->scores[i] = G->deg[in][i] + G->deg[out][i];

  /* s->inCList = (int*)malloc(G->v*sizeof(int)); */
#ifndef NDEBUG
  s->sizeLastBestSketch = 0;
  bitmapAlloc(s->lastVerts, G->v);
#endif

  s->DegEx[in] = allocDA();
  s->DegEx[out] = allocDA();

  s->o = allocO(G, s->t);

  return s;
}
void
freeS(sa_state s)
{
  /* free(s->inCList); */
  freeSTree(s->t);
  freeSTree(s->O[in]);
  freeSTree(s->O[out]);
  free(s->Top2V);
  free(s->condemned);
  free(s->cList2);
  free(s->posInCLs);
  // free(s->cBool);
  free(s->cList);
  free(s->scores);

  freeDA(s->DegEx[in]);
  freeDA(s->DegEx[out]);
  freeO(s->o);

#ifndef NDEBUG
  bitmapDestroy(s->lastVerts);
#endif

  free(s);
}

static void
swapI(int *A,
      int *B
      )
{
  static int T;

  T = *A;
  *A = *B;
  *B = T;
}

static void
shuffle(
  int *a,
  int *idx,
  size_t n
) {
  if (n < 2)
    return; // it is shuffled
  unsigned rnd1;
  for (size_t i = 0; i < n/2; ++i)
  {
    // do
      rnd1 = WELLRNG512a_u4()%(n-i);
    // while (!rnd1);
    idx[a[i]] = i+rnd1;
    idx[a[i+rnd1]] = i;
    swapI(&a[i], &a[i+rnd1]);
  }
  for (size_t i = n/2+1; i < n; ++i)
  {
    // do
      rnd1 = WELLRNG512a_u4()%(n-i);
    // while (!rnd1);
    idx[a[i]] = i-1-rnd1;
    idx[a[i-1-rnd1]] = i;
    swapI(&a[i], &a[i-1-rnd1]);
  }
}

static void
dumpCandidates(sa_state s)
{
  if (s->cLs2 == 0)
    return;
  for (int i = 0; i < s->cLs; ++i)
    swapI(&(s->cList[i]), &(s->cList2[s->cLs2++]));
  for (int i = 0; i < s->cLs2; ++i)
    s->posInCLs[s->cList2[i]] = i;
  int *tmp = s->cList; // swap lists
  s->cList = s->cList2;
  s->cLs = s->cLs2;
  s->cList2 = tmp;
  s->cLs2 = 0; // empty
}

static void
callReorganize(sa_state s)
{
  dumpCandidates(s);
  // call reorganize 
  s->cList[s->cLs] = -1;
  reorganize(s->o, &s->cLs, s->posInCLs, s->cList);
}

void
choose(sa_state s
       )
{
  if (0 >= s->cLs)
    dumpCandidates(s);
  assert(0 < s->cLs && "Empty candidate list");
  // if (0 >= s->cLs) return -1;
  if (s->c < s->cLs)
  {
    s->posInCLs[s->cList[0]] = s->c;
    s->posInCLs[s->cList[s->c]] = 0;
    swapI(&(s->cList[0]), &(s->cList[s->c]));
  }
  else
    s->c = 0;
  s->c++;
}

/* Returns how many elements will be removed if position
   p is used. */
static int
incompatibleSize(sa_state s,
		 int p
		 )
{
  int r = 0;
  node floor;
  node ceil;

  roundSt(s->O[out], p, &floor, &ceil);
  if(NULL != ceil)
    r += value(s->O[out], ceil, NULL);
  else
    r += size(s->O[out]);
  /* Counting strictly less than p */

  roundSt(s->O[in], p, &floor, &ceil);
  if(NULL != floor){
    int c;
    value(s->O[in], floor, &c);
    r += c;

    if(floor == ceil){ /* Also equal to p */
      r++;
    }
  }
  else
    r += size(s->O[in]);
  /* Counting less than or equal to p */

  return r;
}

static int
next(sa_state s,
     int pos,
     int bound
     )
{
  int p = INT_MAX;

  node floor;
  node ceil;

  roundSt(s->O[in], pos, &floor, &ceil);
  if(NULL != ceil)
    p = 1 + key(s->O[in], ceil);

  if(pos < bound && p > bound)
    p = bound;

  return p;
}

static int
prev(sa_state s,
     int pos,
     int bound
     )
{
  int p = -1;
  int pout = -1;

  node floor;
  node ceil;

  roundSt(s->O[out], pos-1, &floor, &ceil);
  if(NULL != floor)
    pout = key(s->O[out], floor);

  if(pout > p)
    p = pout;

  if(pos > bound && p < bound)
    p = bound;

  return p;
}

int
forceReorganize(sa_state s
      )
{
  callReorganize(s);
  // s->c = 0;
}

int
tryReorganize(sa_state s
      )
{
  /* Apply re-organization heuristic */
  if (checkO(s->o, 0))
  {
#if 0
// #ifndef NDEBUG
    int L[s->G->v+1];
    int *R = getInorder(s->t, &L[0], NULL, left);
    *R = -1;
    R = &L[0];
    fprintf(stderr, "  [before-re]maxDAG: ");
    while (-1 != *R)
      fprintf(stderr, "%d ", /*s->G->d[*/*(R++)/*]+1*/);
    fprintf(stderr, "(%ld)\n", R-L);
    assertState(s);
#endif /* NDEBUG */
    callReorganize(s);
#if 0
// #ifndef NDEBUG
    R = getInorder(s->t, &L[0], NULL, left);
    *R = -1;
    R = &L[0];
    fprintf(stderr, "  [after-re] maxDAG: ");
    while (-1 != *R)
      fprintf(stderr, "%d ", /*s->G->d[*/*(R++)/*]+1*/);
    fprintf(stderr, "(%ld)\n", R-L);
    assertState(s);
#endif /* NDEBUG */
  }
}

static void
addIncoming(sa_state s,
      int *pa,
      int *pb,
      int **inA,
      int dE
      )
{
  /* Add one incoming */
  node n;
  if(*pa <= *pb && -1 != **inA){
    if(NULL == (n = getNode(s->t, **inA)))
      pushDA(s->DegEx[in], **inA);
    else { /* Vertex exists in DAG */
      int top = value(s->t, n, NULL);
      s->Top2V[top] = **inA;
      insertKey(s->O[in], top); // stores the order of the node

      /* Move pa forward */
      while(
        *pa <= *pb &&
        incompatibleSize(s, *pa) > 1 - dE
      ) {
        *pa = next(s, *pa, 1+*pb);
      }
    }
    (*inA)++; /* Cycle step */
  }
}

static void
addOutgoing(sa_state s,
      int *pa,
      int *pb,
      int **outA,
      int dE
      )
{
  /* Add one outgoing */
  node n;
  if(*pa <= *pb && -1 != **outA){
    if(NULL == (n = getNode(s->t, **outA)))
      pushDA(s->DegEx[out], **outA);
    else { /* Vertex exists in DAG */
      int top = value(s->t, n, NULL);
      s->Top2V[top] = **outA;
      insertKey(s->O[out], top);

      /* Move pb backward */
      while(
        *pa <= *pb &&
        incompatibleSize(s, *pb) > 1 - dE
      ) {
        *pb = prev(s, *pb, *pa); // count this as processed node
      }
    }
    (*outA)++; /* Cycle step */
  }
}

static void
selectSplitPoint(sa_state s,
      int *pa,
      int *pb
      )
{
  /* Selecting the actual split point */
  s->p = *pa;
  if(*pa <= *pb){ /* Otherwise check failed */
    int count = 0;
    int pmin = INT_MAX;
    int p = *pa;
    while(p <= *pb){
      int incomp = incompatibleSize(s, p);
      if(incomp < pmin){
        pmin = incomp;
        count = 0;
      }
      if(incomp == pmin)
        count++;
      p = next(s, p, 1+*pb);
    }

    if(1 < count){
      count = ARC4RAND_UNIF(count);
      count++;
    }

    p = *pa;
    while(0 < count){
      if(incompatibleSize(s, p) == pmin)
	      count--;
      if(0 < count){
	      p = next(s, p, 1+*pb);
      }
    }

    /* Choose uniformly in interval */
    int n = next(s, p, 1+*pb);
    if (p < *pb && p+1 < n){
      p += ARC4RAND_UNIF(n - p);
    }

    s->p = p;
  }
}

static void
reloadGraph(sa_state s,
      int v
      )
{
  // change the order of the candidates
  int *T = getInorder(s->O[in], s->G->E[v][in], s->Top2V, right);
  dumpDA(s->DegEx[in], T);
#ifndef NDEBUG
  for (int i = 0; i < s->cLs; ++i)
    assert(s->posInCLs[s->cList[i]] == i && "invalid position");
#endif
  /* int c = s->c+1;
  if (s->c & 0xF == 0xF)
  { // reorder candidates not very frequently
    while (-1 != *T)
    { // cannot change the position 0 
      if (c >= s->cLs)
        c = 1;
      int v = *(T++);
      int p = s->posInCLs[v];
      if (0 <= p)
      {
        s->posInCLs[v] = c; 
        s->posInCLs[s->cList[c]] = p; 
        swapI(&(s->cList[c++]), &(s->cList[p]));
      }
    }
  } */
  T = getInorder(s->O[out], s->G->E[v][out], s->Top2V, left);
  dumpDA(s->DegEx[out], T);
  /* if (s->c & 0xF == 0xF)
  {
    while (-1 != *T)
    {
      if (c >= s->cLs)
        c = 1;
      int v = *(T++);
      int p = s->posInCLs[v];
      if (0 <= p)
      {
        s->posInCLs[v] = c; 
        s->posInCLs[s->cList[c]] = p; 
        swapI(&(s->cList[c++]), &(s->cList[p]));
      }
    }
  } */
}

static int
computeCheckDE(sa_state s
      )
{
  return 1 - incompatibleSize(s, s->p);
}

/* Returns the proposed energy variation */
int
check(sa_state s,
      int dE
      )
{
  // printf("cList (size t = %i): ", size(s->t));
  // for (int i = 0; i < s->cLs; ++i) printf("%i ", s->cList[i]);
  // printf("\n");
  assert(2 > dE && "Check with impossible dE");

  /* Candidate */
  int v = s->cList[0];
  int *inA = s->G->E[v][in];
  int *outA = s->G->E[v][out];

  /* printf("State before check\n"); */
  /* printS(s); */
  /* printf(">>>> Candidate is : %d\n", v); */

  /* Load trees O[in] and O[out] */
  clearTree(s->O[in]);
  clearTree(s->O[out]);
  resetDA(s->DegEx[in]);
  resetDA(s->DegEx[out]);

  /* Positions over Topological order */
  int pa = 0;
  int pb = size(s->t);

  while  (pa <= pb &&
        !(-1 == *inA && -1 == *outA)
	){
    /* Increment the heuristic */
    checkO(s->o, 2);

    /* Add one incoming */
    addIncoming(s, &pa, &pb, &inA, dE);

    /* Add one outgoing */
    addOutgoing(s, &pa, &pb, &outA, dE);
  }

  /* Selecting the actual split point */
  selectSplitPoint(s, &pa, &pb);

  /* Reload graph Adj */
  reloadGraph(s, v);

  /* printf("State after check\n"); */
  /* printS(s); */

#ifndef NDEBUG
  assertState(s);
#endif /* NDEBUG */

  s->checkDE = computeCheckDE(s);

  /* fprintf(stderr, "pa = %d ", pa); */
  /* fprintf(stderr, "p = %d ", s->p); */
  /* fprintf(stderr, "pb = %d ", pb); */
  /* fprintf(stderr, "dE = %d ", dE); */
  /* fprintf(stderr, "\n"); */

  // printf("--- checkDE = %i, p = %i\n", s->checkDE, s->p);
  return s->checkDE;
}

/* Executes validated transition */
/* Returns number of new candidates */
void
execute(sa_state s
        )
{
  /* printf("State before execute\n"); */
  /* printS(s); */
  /* fprintf(stderr, "Before candidate insertion"); */
  /* fprintf(stderr, "\n"); */

  int v = s->cList[0]; /* The candidate */
  s->cLs--;
  if (0 < s->cLs)
  {
    s->posInCLs[s->cList[s->cLs]] = 0;
    swapI(&(s->cList[0]), &(s->cList[s->cLs]));
  }
  s->posInCLs[v] = -1;
  /* printf("Adopting candidate %d\n", v); */
#ifdef POSTSCRIPT
  fprintf(LOG, "%d %d i\n", 1+(v%POST_NB_VERT), 1+(v/POST_NB_VERT));
#endif /* POSTSCRIPT */

  // fprintf(stderr, "size_tree = %i, ins %d in pos = %i\n",
  //   size(s->t), v, s->p);

  /* 1 - add v */
  insertInorderKey(s->t, s->p, v);

  assert(size(s->t) > 0);

  /* printS(s); */
  /* fprintf(stderr, "After candidate insertion"); */
  /* fprintf(stderr, "\n"); */

  // printf("checkDE = %i, incompatible size = %i, p = %i, size_tree = %i\n",
  //   s->checkDE, incompatibleSize(s, s->p), s->p, size(s->t)/*, s->t->root*/);

  /* 2 - remove incompatible nodes */
  int *condemned = s->condemned; /* Relabel and re-use */
  splitSt(s->O[in], s->p, right);
  // TODO: is this needed? there are repeated nodes in condemned...
  condemned = getInorder(s->O[in], condemned, s->Top2V, right);

  splitSt(s->O[out], s->p-1, left);
  condemned = getInorder(s->O[out], condemned, s->Top2V, left);
  *condemned = -1;

  condemned = s->condemned;
  while(-1 != *condemned && 0 >= s->checkDE){
    removeN(s->t, getNode(s->t, *condemned));
    /* Return back to candidate list, as current candidate. */
    /* printf("Returning candidate %d\n", s->Top2V[*condemned]); */
#ifdef POSTSCRIPT
    fprintf(LOG, "%d %d c\n",
	    1+(*condemned%POST_NB_VERT),
	    1+(*condemned/POST_NB_VERT));
#endif /* POSTSCRIPT */

    /* printS(s); */
    // fprintf(stderr, "Removed node %d\n", *condemned);

    s->cList2[s->cLs2] = *condemned;
    s->cLs2++;
    condemned++;
    s->checkDE++;
  }

#ifndef NDEBUG
  assertState(s);
#endif /* NDEBUG */
}

void
prepareS_vertCover(sa_state s,
         int*  A, /* Array with vertexes */
         int   l	/* length of the array */
) {
  s->cLs = 0;
  clearTree(s->t);
  s->c = 0;

  /* Dump array into the tree */
  for (int i = 0; i < l; i++)
    reRoot(s->t, A[i], left);

  /* Run approximation discarding */
  for(int i = 0; i < l; i++){
    int v = A[i];
    /* Identify position in topological order */
    node nv = getNode(s->t, v);
    if(NULL != nv){
      int top = value(s->t, nv, NULL);
      int *outA = s->G->E[v][out];

      node u;
      int valid = 1;
      while (-1 != *outA && valid) {
        // s->cBool[*outA] = 1;
        u = getNode(s->t, *outA);
        if (NULL != u &&
            value(s->t, u, NULL) <= top
        ) valid = 0;
        outA++;
        // nbProcVert++;
      }

      if (!valid) { /* Means early temination */
        removeN(s->t, u);
        s->cList[s->cLs] = outA[-1];
        s->cLs++;
        removeN(s->t, nv);
        s->cList[s->cLs] = v;
        s->cLs++;
      }
    }
  }
  shuffleCandidates(s);
}

int
getNbCandidates(sa_state s
      )
{
  return s->cLs;
}

void
tryAllCandidates(sa_state s,
  int pdE
      )
{
  int r = 0;
  for (int j = 0; j < s->cLs; j++)
  {
    if (1 == check(s, pdE))
      execute(s);
    ++r; // put last in r-th place (r is the next in the initial order)
    if (r < s->cLs)
    {
      s->posInCLs[s->cList[0]] = r;
      s->posInCLs[s->cList[r]] = 0;
      swapI(&(s->cList[0]), &(s->cList[r]));
    }
  }
  dumpCandidates(s);
}

void
shuffleCandidates(sa_state s
      )
{
  s->c = 0;
  assert(0 != s->cLs && "empty candidates");
  shuffle(s->cList, s->posInCLs, s->cLs);
}

void
decayAndReorganize(sa_state s
      )
{
  /* tryAllCandidates(s, 0); */
  for (int i = 0; i < 1; ++i)
  {
    s->cList[s->cLs] = -1;
    callReorganize(s);
    tryAllCandidates(s, 1);
  }
}

void
prepareS_greedy(sa_state s,
         int*  A, /* Array with vertexes */
         int   l	/* length of the array */
) {
  s->cLs = 0;
  clearTree(s->t);
  s->c = 0;

  // for (int i = l-1; 0 <= i; i--) {
  for (int i = 0; i < l; i++) {
    s->cList[s->cLs] = A[i];
    s->posInCLs[A[i]] = s->cLs;
    s->cLs++;
  }
  // sortVertexArrayOnDegProd_largeLast(s->G, s->cList, s->cLs);
  // qsort_r(s->cList, s->cLs, sizeof(int), pcmp, s->G);
  
  tryAllCandidates(s, 0);
  callReorganize(s);
  tryAllCandidates(s, 1);
}

void
prepareS_empty(sa_state s
      ) {
  s->cLs = 0;
  clearTree(s->t);
  s->c = 0;

  for (int i = 0; i < s->G->v; i++)
  {
    s->posInCLs[i] = s->cLs;
    s->cList[s->cLs++] = i;
  }
  // sortVertexArrayOnDegProd_largeFirst(s->G, s->cList, s->cLs);
  // qsort_r(s->cList, s->cLs, sizeof(int), pcmp, s->G);
  shuffleCandidates(s);
  callReorganize(s);
}

void
prepareS(sa_state s,
         int*  A, /* Array with vertexes */
         int   l	/* length of the array */
) {
  s->cLs = 0;
  clearTree(s->t);
  s->c = 0;

#ifndef NUSE_INIT_APPROX
  /* Dump array into the tree */
  for (int i = 0; i < l; i++) {
    reRoot(s->t, A[i], left);
    // s->cBool[A[i]] = 1;

#ifdef POSTSCRIPT
    fprintf(LOG, "%d %d i\n",
	    1+(A[i]%POST_NB_VERT),
	    1+(A[i]/POST_NB_VERT));
    fflush(LOG);
#endif /* POSTSCRIPT */
  }

  /* Run approximation discarding */
  for(int i = 0; i < l; i++){
    int v = A[i];
    /* Identify position in topological order */
    node nv = getNode(s->t, v);
    if(NULL != nv){
      int top = value(s->t, nv, NULL);
      int *outA = s->G->E[v][out];

      node u;
      int valid = 1;
      while (-1 != *outA && valid) {
        // s->cBool[*outA] = 1;
        u = getNode(s->t, *outA);
        if (NULL != u &&
            value(s->t, u, NULL) <= top
        ) valid = 0;
        outA++;
        // nbProcVert++;
      }

      if (!valid) { /* Means early temination */
        removeN(s->t, u);
        s->cList[s->cLs] = outA[-1];
        s->cLs++;
#ifdef POSTSCRIPT
	fprintf(LOG, "%d %d c\n",
		1+(outA[-1]%POST_NB_VERT),
		1+(outA[-1]/POST_NB_VERT)); // TODO
	fflush(LOG);
#endif /* POSTSCRIPT */
        removeN(s->t, nv);
        s->cList[s->cLs] = v;
        s->cLs++;
#ifdef POSTSCRIPT
	fprintf(LOG, "%d %d c\n",
		1+(v%POST_NB_VERT),
		1+(v/POST_NB_VERT));
	fflush(LOG);
#endif /* POSTSCRIPT */
      }
    }
  }
#else
  for (int i = 0; i < l; i++) {
    s->cList[s->cLs] = A[i];
    s->cLs++;
  }
#endif /* NUSE_INIT_APPROX */

  dumpCandidates(s);
  shuffleCandidates(s);
}

int
getE(sa_state s)
{
  return size(s->t);
}

#ifndef NDEBUG
static void
gdbBreak(void)
{}

static void
assertState(sa_state s)
{
// #if 0
  for (int i = 0; i < s->G->v; i++)
  {
    node u = getNode(s->t, i);
    if (NULL != u)
    {
      int posU = value(s->t, u, NULL);
      int *outA = s->G->E[i][out];
      while (-1 != *outA)
      {
        node v = getNode(s->t, *outA);
        if (NULL != v)
        {
          int posV = value(s->t, v, NULL);
          // fprintf(stderr, "test edge %i(pos=%i)->%i(pos=%i) \n", i, posU, *outA, posV);
          if(posU >= posV)
            gdbBreak();
          assert(posU < posV && "Invalid sa_state configuration");
        }
        outA++;
      }
    }
  }
// #endif
}
#endif /* NDEBUG */

/* Modified sa_state print for torus */

void
printS(sa_state s)
{
  int side = ceil(sqrt(s->G->v));
  int *A = malloc(side*side*sizeof(int));

  for(int i = 0; i < side*side; i++)
    A[i] = -1;

  if(NULL != s->t){
    int *L = (int *)malloc((1+size(s->t))*sizeof(int));
    *getInorder(s->t, &L[0], NULL, left) = -1;

    int j = 0;
    int *T = L;
    while(-1 != *T){
      A[*T] = j++;
      T++;
    }
    free(L);
  }

  int v = s->cList[0];
  for(int i = 0; i < side*side; i++){
    if(v == i) {
      fprintf(stderr, " * ");
    } else if(-1 != A[i]) {
      fprintf(stderr, "%.2d ", A[i]);
    } else {
      fprintf(stderr, "   ");
    }
    if(0 == ((1+i)%side)){
      fprintf(stderr, "\n");
    }
  }

  free(A);
}

void
printS_v2(sa_state s)
{
  // printf("Printing sa_state\n");
  if(NULL != s->t){
    int *L = (int *)malloc((1+size(s->t))*sizeof(int));
    int *T = getInorder(s->t, L, NULL, left);
    *T = -1;

    fprintf(stderr, "Order (size=%ld): ", T-L);
    int *t = L;
    while(-1 != *t && 3 > t-L){
      fprintf(stderr, "%d ", *t);
      t++;
    }
    fprintf(stderr, "... ");
    fprintf(stderr, "%d\n", *(T-1));
    free(L);
  }

  fprintf(stderr, "1st candidate: %i, last candidate: %i\n", s->cList[0], s->cList[s->cLs-1]);
}

int *
getSketch(sa_state s,
	        int *L
          )
{
  int *R = getInorder(s->t, L, NULL, left);
  *R = -1;

#if 0
// #ifndef NDEBUG
  int sizeSketch = R - L;
  int *l = L;
  int diff = 0;
  bitmap newVerts;
  bitmapAlloc(newVerts, s->G->v);
  fprintf(stderr, "getSketch diff: ");
  if (sizeSketch > s->sizeLastBestSketch)
  {
    while (-1 != *l)
    {
      if (!bitmapCheck(s->lastVerts, *l))
      {
        fprintf(stderr, "%d[%d] ", *l, s->G->deg[in][*l]*s->G->deg[out][*l]);
        diff++;
      }
      bitmapPut(newVerts, *l);
      l++;
    }
  }
  fprintf(stderr, "(sOld=%d sNew=%d diff=%d)\n", s->sizeLastBestSketch, sizeSketch, diff);
  s->sizeLastBestSketch = sizeSketch;
  bitmapDestroy(s->lastVerts);
  s->lastVerts = newVerts;
#endif

  return R;
}
