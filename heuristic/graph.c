#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-value"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"

#include "bitmap.h"
#include "graph.h"
#define PRINT_GRAPH 0

#if PRINT_GRAPH
static int alreadyPrinted = 0;
#endif

static void
tarjanAllocG(graph G, int inPlace);

static void
tarjanFreeG(graph G);

graph
allocG(long v, long e)
{
  graph G = (graph)calloc(1, sizeof(struct graph_));
  G->v = v;
  G->e = e;
  G->deg[in]  = calloc(v, sizeof(int));
  G->deg[out] = calloc(v, sizeof(int));
  G->E = malloc(v*2*sizeof(int*));
  G->t = NULL;
  return G;
}

/* static int
sortById(
  const void *p1, 
  const void *p2
) {
  int i1 = *(int *)p1;
  int i2 = *(int *)p2;
  return i2 - i1;
} */

static void 
giveScores(
  graph G,
  int v,
  char *L,
  int *score
) {
  if (!L[v])
    L[v] = 1;
  int *o = G->E[v][out];
  while (-1 != *o)
  {
    if (!L[v])
      giveScores(G, *o, L, score);
    else
      score[*o]++; // the intersection vertex gets the score
    o++;
  }
}

void
scoreVertexByLoopSize(
  graph G
) {
  char *L = (char*)malloc(G->v*sizeof(char));
  if (!G->score)
    G->score = (int*)calloc(G->v, sizeof(int));
  for (int v = 0; v < G->v; v++)
  {
    memset(L, 0, G->v*sizeof(char));
    if (!G->score[v])
      giveScores(G, v, L, G->score);
  }

  free(L);
}

graph
extractSCC(
  graph G,
  int *sccList,
  int lenScc
) {
  int *p = sccList, *ptr, *o, *i;
  int it = 0;
  int *translate, *d;
  int *invTranslate = G->d;
  int countE = 0;
  graph G1;
  translate = d = (int*)malloc(lenScc*sizeof(int));
#ifndef NDEBUG
  memset(invTranslate, -1, lenScc*sizeof(int));
#endif
  while (it < lenScc)
  {
    invTranslate[*p] = it;
    *d = *p;
    // printf("%i -> %i\n", it, *p);
    assert(-1 != invTranslate[*p] && "invalid vertex in SCC");
    o = G->E[*p][out];
    while (-1 != *o) { countE++; o++; }
    d++; p++; it++;
  }
  G1 = allocG(lenScc, countE);
  G1->d = translate;
  ptr = malloc(2*(lenScc+countE+1)*sizeof(int));
  p = sccList;
  it = 0;
  while (it < lenScc)
  {
    /* handle out vertices */
    o = G->E[*p][out];
    G1->E[invTranslate[*p]][out] = ptr;
    // int count = 0;
    while (-1 != *o)
    {
      assert(-1 != invTranslate[*o] && "invalid vertex in SCC adjacency list");
      *ptr = invTranslate[*o]; ptr++; o++;
      // count++;
      G1->deg[out][invTranslate[*p]]++;
    }
    // qsort(G1->E[invTranslate[*p]][out], count, sizeof(int), sortById); // TODO
    *ptr = -1; ptr++;

    /* handle in vertices */
    i = G->E[*p][in];
    G1->E[invTranslate[*p]][in] = ptr;
    // count = 0;
    while (-1 != *i)
    {
      assert(-1 != invTranslate[*i] && "invalid vertex in SCC adjacency list");
      *ptr = invTranslate[*i]; ptr++; i++;
      G1->deg[in][invTranslate[*p]]++;
      // count++;
    }
    // qsort(G1->E[invTranslate[*p]][in], count, sizeof(int), sortById); // TODO
    *ptr = -1; ptr++;
    p++; it++;
  }

  for (int v = 0; v < G1->v; ++v)
  {
    sortVertexArrayOnDegProd_largeFirst(G1, G1->E[v][out], G1->deg[out][v]);
    sortVertexArrayOnDegProd_largeFirst(G1, G1->E[v][in], G1->deg[in][v]);
  }

  return G1;
}

void
freeG(graph G)
{
  if (G->E)
  {
    free(G->start_E);
    free(G->E);
    G->E = NULL;
  } 
  if (G->d) free(G->d);
  if (G->Eaux) 
  {
    free(G->start_Eaux);
    free(G->Eaux);
    G->E = NULL;
  }
  if (G->inFVS)
    free(G->inFVS);
  if (G->score) free(G->score);
  tarjanFreeG(G);
  free(G);
}

// matrix accessed via idx_row*nbVertices+idx_col
graph
fromSquareMat(long nbVertices, unsigned char *mat)
{
  int edgeEstimate = nbVertices;
  graph G;

  int (*E)[2] = (int (*)[2])malloc(edgeEstimate*2*sizeof(int*));
  int countE = 0;
  for(int i = 0; i < nbVertices; i++){
    for(int j = 0; j < nbVertices; j++){
      int idx = i * nbVertices + j;
      if (mat[idx] != 0) {
        E[countE][out] = i;
        E[countE][in]  = j;
        assert(E[countE][out] != E[countE][in] && "Loop node in input");
        countE++;
        if (countE > edgeEstimate) {
          edgeEstimate <<= 1;
          E = realloc(E, edgeEstimate*2*sizeof(int*));
        }
      }
    }
  }
  G = allocG(nbVertices, countE);

  /* Count Sort for packing */
  int *c = calloc(2*G->v, sizeof(int)); /* Clean-up */
  for(int i=0; i<G->e; i++){ /* Count */
    int outE = 2*E[i][out];
    int inE = 2*E[i][in];
    c[out+outE]++;
    c[in+inE]++;
  }

  /* Pointer */
  int *p = malloc(2*(G->v+G->e)*sizeof(int));
  for(int i=0; i<2*G->v; i++){ /* Acc */
    p += c[i];
    if(out == i%2)
      G->E[i/2][out] = p;
    else
      G->E[i/2][in] = p;
    *p = -1; /* Add terminator */
    p++;
  }
  free(c);

  G->deg[out] = (int*)calloc(G->v, sizeof(int));
  G->deg[in] = (int*)calloc(G->v, sizeof(int));
  for(int i=0; i<G->e; i++){ /* Transfer */
    G->E[E[i][out]][out]--;
    *(G->E[E[i][out]][out])=E[i][in];
    G->deg[out][E[i][out]]++;
    G->E[E[i][in]][in]--;
    *(G->E[E[i][in]][in])=E[i][out];
    G->deg[in][E[i][in]]++;
  }
#ifndef NDEBUG
  for (int i = 0; i < G->v; i++) {
    assert(G->deg[in][i] == degree(G, i, in));
    assert(G->deg[out][i] == degree(G, i, out));
  }
#endif
  free(E);

  return G;
}

graph
loadG(FILE *stream)
{
  graph G = (graph)calloc(1, sizeof(struct graph_));
  int weights, i, countEdges = 0, countVerts = 0;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;

  while ((read = getline(&line, &len, stream)) != -1)
    if (line[0] != '%')
      break;
  
  sscanf(line, "%ld %ld %d\n", &(G->v), &(G->e), &weights);
  // assert(weights == 0);
  
  G->t = NULL;

  G->E = malloc(G->v*2*sizeof(int *));

  int (*E)[2] = (int (*)[2])malloc(G->e*2*sizeof(int));

  while ((read = getline(&line, &len, stream)) != -1) {
    if (line[0] == '%') continue; // ignore
    if (line[0] == '\n') 
    {
      countVerts++;
      continue; // ignore
    }
    // fprintf(stderr, "before strtok 1\n");
    char *vert = strtok(line, " ");
    // fprintf(stderr, "after strtok 1\n");
    int adjV;
    while (vert != NULL && vert[0] != '\n') {
      assert(countEdges < G->e && "The number of edges does not match input");
      // fprintf(stderr, "before atoi 1\n");
      adjV = atoi(vert);
      // fprintf(stderr, "after atoi 1\n");
      E[countEdges][out] = countVerts;
      E[countEdges][in] = adjV-1;
      assert(E[countEdges][out] != E[countEdges][in] && "Loop node in input");
      assert(E[countEdges][out] >= 0 && E[countEdges][out] < G->v && "Invalid vertex ID");
      assert(E[countEdges][ in] >= 0 && E[countEdges][ in] < G->v && "Invalid vertex ID");
      // fprintf(stderr, "before strtok 3\n");
      vert = strtok(NULL, " ");
      // fprintf(stderr, "after strtok 3\n");
      countEdges++;
    }
    countVerts++;
  }
  assert(countEdges == G->e && "The number of edges does not match input");
  assert(countVerts == G->v && "The number of vertices does not match input");

  /* Count Sort for packing */
  int *c = calloc(2*G->v, sizeof(int)); /* Clean-up */
  for(int i=0; i<G->e; i++){ /* Count */
    c[out+2*E[i][out]]++;
    c[in+2*E[i][in]]++;
  }

  /* Pointer */
  int *p = malloc(2*(G->v+G->e)*sizeof(int));
  for(i=0; i<2*G->v; i++){ /* Acc */
    p += c[i];
    if(out == i%2)
      G->E[i/2][out] = p;
    else
      G->E[i/2][in] = p;
    *p = -1; /* Add terminator */
    p++;
  }
  free(c);

  G->deg[out] = (int*)calloc(G->v, sizeof(int));
  G->deg[in] = (int*)calloc(G->v, sizeof(int));
  for(i=0; i<G->e; i++){ /* Transfer */
    G->E[E[i][out]][out]--;
    *(G->E[E[i][out]][out])=E[i][in];
    G->deg[out][E[i][out]]++;
    G->E[E[i][in]][in]--;
    *(G->E[E[i][in]][in])=E[i][out];
    G->deg[in][E[i][in]]++;
  }
#ifndef NDEBUG
  for (i = 0; i < G->v; i++) {
    assert(G->deg[in][i] == degree(G, i, in));
    assert(G->deg[out][i] == degree(G, i, out));
  }
#endif
  free(E);

  return G;
}

void
setEaux(graph G)
{
  G->Eaux = (int *(*)[2])malloc(2*G->v*sizeof(int*));
  int *p = (int*)malloc((G->v+G->e)*2*sizeof(int));
  G->start_Eaux = p;

  intptr_t diff_base = ((intptr_t)G->start_Eaux) - ((intptr_t)G->start_E);
  for (int v = 0; v < G->v; ++v)
  {
    intptr_t po_orig = (intptr_t)(G->E[v][out]);
    intptr_t pi_orig = (intptr_t)(G->E[v][in]);
    po_orig += diff_base;
    pi_orig += diff_base;
    G->Eaux[v][out] = (int*)po_orig;
    G->Eaux[v][in] = (int*)pi_orig;
    int *i = G->E[v][in], *iaux = G->Eaux[v][in], *o = G->E[v][out], *oaux = G->Eaux[v][out];
    while (*(i++) != -1)
      *(iaux++) = 0;
    *iaux = -1;
    while (*(o++) != -1)
      *(oaux++) = 0;
    *oaux = -1;
  }
}

#if PRINT_GRAPH
graph
loadG(FILE *stream)
{
  graph G = (graph) malloc(sizeof(struct graph_));
  fscanf(stream, "%ld%ld", &(G->v), &(G->e));
  G->t = NULL;

// #ifndef NDEBUG
  if (!alreadyPrinted) {
    printf("# V=%li,E=%li,d(v)=%f,logV=%f\n", G->v, G->e, (float)G->e/(float)G->v, log2f(G->v));
    alreadyPrinted = 1;
  }
// #endif

  G->E = malloc(G->v*2*sizeof(int *));

  int (*E)[2] = (int (*)[2])malloc(G->e*2*sizeof(int));
  for(int i = 0; i < G->e; i++){
    // TODO: check input issues like repeated edges (self-loops fails the assert)
    fscanf(stream, "%d%d", &(E[i][out]), &(E[i][in]));

    // TODO: add these 2 lines for TangEtAl2017 graphs (vertex ID starts at 1)
    // (E[i][out])--;
    // (E[i][in])--;

    assert(E[i][out] != E[i][in] && "Loop node in input");
    assert(E[i][out] >= 0 && E[i][out] < G->v && "Invalid vertex ID");
    assert(E[i][ in] >= 0 && E[i][ in] < G->v && "Invalid vertex ID");
  }

  /* Count Sort for packing */
  int *c = calloc(2*G->v, sizeof(int)); /* Clean-up */
  for(int i=0; i<G->e; i++){ /* Count */
    c[out+2*E[i][out]]++;
    c[in+2*E[i][in]]++;
  }

  /* Pointer */
  int *p = malloc(2*(G->v+G->e)*sizeof(int)); // TODO: shouldn't be (2*G->v+G->e) ?
  for(int i=0; i<2*G->v; i++){ /* Acc */
    p += c[i];
    if(out == i%2)
      G->E[i/2][out] = p;
    else
      G->E[i/2][in] = p;
    *p = -1; /* Add terminator */
    p++;
  }
  free(c);

  G->deg[out] = (int*)calloc(G->v, sizeof(int));
  G->deg[in] = (int*)calloc(G->v, sizeof(int));
  for(int i=0; i<G->e; i++){ /* Transfer */
    G->E[E[i][out]][out]--;
    *(G->E[E[i][out]][out])=E[i][in];
    G->deg[out][E[i][out]]++;
    G->E[E[i][in]][in]--;
    *(G->E[E[i][in]][in])=E[i][out];
    G->deg[in][E[i][in]]++;
  }
#ifndef NDEBUG
  for (int i = 0; i < G->v; i++) {
    assert(G->deg[in][i] == degree(G, i, in));
    assert(G->deg[out][i] == degree(G, i, out));
  }
#endif
  free(E);

  return G;
}
#endif

void
printG(graph G)
{
  for(int i = 0; i < G->v; i++){
    printf("vertex %d\n", i);
    printf("out: ");
    int *p = G->E[i][out];
    while(-1 != *p){
      printf("%d ", *p);
      p++;
    }
    printf("\n");

    printf("in: ");
    p = G->E[i][in];
    while(-1 != *p){
      printf("%d ", *p);
      p++;
    }
    printf("\n");
  }
}

static int countEdgesPerSCC(graph G, int count_rG)
{
  int *v = &(G->t->order[1]);
  int *scc = G->t->SCCsizeS;
  int *edges = G->t->SCCnbEdges;
#ifndef NDEBUG
  int count_scc_size = 0;
  while (*scc != -1) {
    count_scc_size += *scc;
    scc++;
  }
  assert(count_scc_size == G->v);
  scc = G->t->SCCsizeS;
#endif
  if (edges) { // optional: count edges per SCC
    while (*scc != -1) {
      *edges = 0;
      for(int i = 0; i < *scc; ++i, ++v) {
        int *outE;
        if (count_rG)
          outE = G->t->rG->E[*v][in];
        else
          outE = G->E[*v][in];
        while(*outE != -1) {
          (*edges)++;
          outE++;
        }
      }
      edges++;
      scc++;
    }
    *edges = -1;
  } else {
    return -1;
  }
  return 0;
}

// Needs a list to store temporary inputs and outputs
typedef struct contract_tmp_list_node_
{
  int v;
  struct contract_tmp_list_node_ *next;
} *
cTmpListN;

typedef struct contract_tmp_list_
{
  int nbV;
  cTmpListN *V[2];
  int *sizes[2];
} *
cTmpList;

static void
cTmpListAlloc(cTmpList *lst, int v)
{
  *lst = calloc(1, sizeof(struct contract_tmp_list_));
  (*lst)->nbV = v;
  (*lst)->V[out] = (cTmpListN*)calloc(v, sizeof(cTmpListN));
  (*lst)->V[in] = (cTmpListN*)calloc(v, sizeof(cTmpListN));
  (*lst)->sizes[out] = (int*)calloc(v, sizeof(int));
  (*lst)->sizes[in] = (int*)calloc(v, sizeof(int));
}

static void
cTmpListFree_rec(cTmpListN n)
{
  if (n->next)
    cTmpListFree_rec(n->next);
  free(n);
}

static void
cTmpListFree(cTmpList lst, int afterFlush)
{
  if (!afterFlush)
  {
    for (int v = 0; v < lst->nbV; ++v)
    {
      if (lst->V[in][v])
        cTmpListFree_rec(lst->V[in][v]);
      if (lst->V[out][v])
        cTmpListFree_rec(lst->V[out][v]);
    }
  }
#ifndef NDEBUG
  else
    for (int v = 0; v < lst->nbV; ++v)
    {
      assert(!lst->V[out][v] && "list is not empty");
      assert(!lst->V[in][v] && "list is not empty");
    }
#endif
  free(lst->V[in]);
  free(lst->V[out]);
  free(lst->sizes[in]);
  free(lst->sizes[out]);
  free(lst);
}

static void
cTmpList_add(cTmpList lst,
  enum adjacencies dir,
  int u,
  int v
      )
{
  assert(0<=v&&v<lst->nbV && "invalid vert");
  assert(0<=u&&u<lst->nbV && "invalid vert");
  cTmpListN n = lst->V[dir][u];
  cTmpListN n2 = (cTmpListN)malloc(sizeof(struct contract_tmp_list_node_));

  n2->v = v;
  n2->next = n;
  lst->V[dir][u] = n2;
  lst->sizes[dir][u]++;
}

static void
cTmpList_checkAndAdd(cTmpList lst,
  enum adjacencies dir,
  int u,
  int v
      )
{
  assert(0<=v&&v<lst->nbV && "invalid vert");
  assert(0<=u&&u<lst->nbV && "invalid vert");
  cTmpListN n = lst->V[dir][u];
  cTmpListN prev = n;

  if (!n)
    return cTmpList_add(lst, dir, u, v);

  while (n)
  {
    if (n->v == v)
      return;
    else
      prev = n;
    n = n->next;
  }

  assert(prev && "cannot insert in this position");
  lst->sizes[dir][u]++;
  prev->next = (cTmpListN)malloc(sizeof(struct contract_tmp_list_node_));
  prev->next->v = v;
  prev->next->next = n;
}

static void
cTmpList_replace(cTmpList lst,
  enum adjacencies dir,
  int u,
  int v,
  int w
      )
{
  assert(0<=v&&v<lst->nbV && "invalid vert");
  assert(0<=u&&u<lst->nbV && "invalid vert");
  assert(0<=w&&w<lst->nbV && "invalid vert");
  cTmpListN n = lst->V[dir][u];

  while (n)
  {
    if (n->v == w)
      break;
    if (n->v == v)
    {
      // fprintf(stderr, "[%d(%s)] replaced %d with %d\n", u, dir==out?"out":"in", v, w);
      n->v = w;
      break;
    }
    n = n->next;
  }

  if (n)
  {
    cTmpListN prev = n;
    n = n->next;
    while (n)
    {
      if (n->v == v || n->v == w)
      { // repeated
        // fprintf(stderr, "[%d(%s)] parallel edge: %d\n", u, dir==out?"out":"in", n->v);
        cTmpListN rm = n;
        prev->next = n->next;
        lst->sizes[dir][u]--;
        free(rm);
        n = prev;
      }
      else
        prev = n;
      n = n->next;
    }
  }
}

static void
cTmpList_remove(cTmpList lst,
  enum adjacencies dir,
  int u,
  int v
      )
{
  assert(0<=u&&u<lst->nbV && "invalid vert");
  assert(0<=v&&v<lst->nbV && "invalid vert");
  cTmpListN n = lst->V[dir][u];
  cTmpListN prev = n;

  if (n && n->v == v)
  {
    // fprintf(stderr, "[%d(%s)] removed %d\n", u, dir==out?"out":"in", v);
    lst->V[dir][u] = n->next;
    free(n);
    lst->sizes[dir][u]--;
    return;
  }

  while (n)
  {
    if (n->v == v)
    {
      // fprintf(stderr, "[%d(%s)] removed %d\n", u, dir==out?"out":"in", v);
      prev->next = n->next;
      free(n);
      lst->sizes[dir][u]--;
      return;
    }
    else
      prev = n;
    n = n->next;
  }
}

static void
cTmpList_removeAllAdj(cTmpList lst,
  enum adjacencies dir,
  int u
      )
{
  assert(0<=u&&u<lst->nbV && "invalid vert");
  cTmpListN n = lst->V[dir][u];
  cTmpListN prev = n;
  while (n)
  {
    cTmpList_remove(lst, OTHER_ADJ(dir), n->v, u);
    n = n->next;
    free(prev);
    prev = n;
  }
  lst->sizes[dir][u] = 0;
  lst->V[dir][u] = NULL;
}

#ifndef NDEBUG
static int countEdgesIn = 0;
#endif

static int *
cTmpList_flushIntoAdjArray(cTmpList lst,
  graph G,
  enum adjacencies dir,
  int u,
  int *Adj
      )
{
  assert(0<=u&&u<lst->nbV && "invalid vert");
  cTmpListN n = lst->V[dir][u];
  cTmpListN prev = n;
  G->deg[dir][u] = lst->sizes[dir][u];
  G->E[u][dir] = Adj;
  if (out == dir)
    G->e += G->deg[dir][u];
#ifndef NDEBUG
  else
    countEdgesIn += G->deg[dir][u];
#endif

  while (n)
  {
    int w = n->v;
    *(Adj++) = w;
    n = n->next;
    free(prev);
    prev = n;
  }
  *(Adj++) = -1;
  lst->V[dir][u] = NULL;
  lst->sizes[dir][u] = 0;
  return Adj;
}

static void
contract1(graph G,
  cTmpList lst,
  enum adjacencies dir,
  int u
      )
{
  assert(1 == lst->sizes[dir][u] && "Invalid out degree");
  assert(NULL == lst->V[dir][u]->next && "Invalid out adj list");
  int v = lst->V[dir][u]->v;
  cTmpListN n = lst->V[OTHER_ADJ(dir)][u];
  cTmpListN nV = lst->V[dir][v];
  // fprintf(stderr, "contracted %d with %d\n", u, v);

  // if out_edges of v points to u, there is a self-loop
  while (nV)
  {
    if (nV->v == u)
    {
      G->inFVS[G->inFVSs++] = v;
      // remove v
      // fprintf(stderr, "put %d in FVS\n", v);
      cTmpList_removeAllAdj(lst, out, v);
      cTmpList_removeAllAdj(lst, in, v);
      goto END; // also remove u
    }
    nV = nV->next;
  }

  while (n)
  {
    // add v in in_edges of u
    cTmpList_checkAndAdd(lst, OTHER_ADJ(dir), v, n->v);
    // replace u in out_edges with v
    cTmpList_replace(lst, dir, n->v, u, v);
    n = n->next;
  }

END:
  // remove u
  cTmpList_removeAllAdj(lst, out, u);
  cTmpList_removeAllAdj(lst, in, u);
}

int
contractG(graph G
      )
{
  int rep = 0;
  int oldFVSs = G->inFVSs;
  static cTmpList lst = NULL;
  int *p = G->start_E;
  if (!G->inFVS)
    G->inFVS = (int*)malloc((G->v+1)*sizeof(int));
  
  if (!lst) // TODO: never freed
    cTmpListAlloc(&lst, G->v);

  for (int v = 0; v < G->v; ++v)
  { // dump vertices into linked list structure
    for(enum adjacencies dir = out; dir <= in; dir++)
    {
      int *i = G->E[v][dir];
      while (-1 != *i)
        cTmpList_add(lst, dir, v, *(i++));
    }
  }

  for (int v = 0; v < G->v; ++v) // contractions
    for(enum adjacencies dir = out; dir <= in; dir++)
      if (1 == lst->sizes[dir][v])
        contract1(G, lst, dir, v);

  G->e = 0;
#ifndef NDEBUG
  countEdgesIn = 0;
#endif
  for (int v = 0; v < G->v; ++v) // merge in the vert
    for(enum adjacencies dir = out; dir <= in; dir++)
      p = cTmpList_flushIntoAdjArray(lst, G, dir, v, p);

  assert(countEdgesIn == G->e && "wrong count of edges");

  if (oldFVSs < G->inFVSs)
    rep = 1;
  return rep;
}

static void
tarjanAllocG(graph G, int inPlace)
{
  if (!G->t) {
    G->t = (TarjanMetadata) calloc(1, sizeof(struct TarjanMetadata_));
    if (!inPlace) G->t->rG = allocG(G->v, G->e);
    else          G->t->rG = NULL;
    G->t->TSbool = (char *)calloc(G->v, sizeof(char));
    G->t->d = (int *)calloc(G->v, sizeof(int));
    G->t->TS = (int *)malloc((1 + G->v + G->e)*sizeof(int)); // TODO: overflow!!!
    G->t->low = (int *)malloc(G->v*sizeof(int));
    G->t->L = (int *)malloc((1+G->v)*sizeof(int));
    *G->t->L = -1;
    // for (int i = 0; i < G->v; ++i) G->t->L[i+1] = i;
    // in calloc
    /* G->t->tSize = 0;
    G->t->visited = 0;
    G->t->top = 0;
    G->t->t = NULL; */
  }
}

static void
tarjanReset(graph G)
{
  memset(G->t->TSbool, 0, G->v*sizeof(char));
  memset(G->t->d, 0, G->v*sizeof(int));
  G->t->top = 0;
  G->t->visited = 0;
}

static void
tarjanFreeG(graph G)
{
  if (G->t) {
    if (G->t->rG)
    {
      // G->E is shared
      if (G->t->rG->d) free(G->t->rG->d);
      if (G->t->rG->Eaux) 
      {
        free(G->t->rG->start_Eaux);
        free(G->t->rG->Eaux);
        G->t->rG->E = NULL;
      }
      if (G->t->rG->score) free(G->t->rG->score);
      if (G->t->rG->inFVS)
        free(G->t->rG->inFVS);
      free(G->t->rG);
    }
    free(G->t->TSbool);
    free(G->t->TS);
    free(G->t->low);
    free(G->t->d);
    free(G->t->L);
    if (G->t->t) 
      free(G->t->t);
    free(G->t);
  }
}

static graph tmp_G;
static void
TarjanVisit(int u
	    )
{
  tmp_G->t->visited++;
  tmp_G->t->d[u] = tmp_G->t->visited;
  tmp_G->t->low[u] = tmp_G->t->visited;

  /* Add discovery items */
  tmp_G->t->TS[tmp_G->t->top] = u;
  tmp_G->t->top++; /* Push u to TS */
  assert(tmp_G->t->top < 1+tmp_G->v+tmp_G->e && "tmp_G->t->top overflows!");
  tmp_G->t->TSbool[u] = 1;

  int *pv = tmp_G->E[u][out];
  while(-1 != *pv){
    if(0 == tmp_G->t->d[*pv])
      TarjanVisit(*pv);
    if(1 == tmp_G->t->TSbool[*pv] && tmp_G->t->low[*pv] < tmp_G->t->low[u])
      tmp_G->t->low[u] = tmp_G->t->low[*pv];

    pv++;
  }

  /* Add finishing items */
  tmp_G->t->TS[tmp_G->t->top] = u;
  tmp_G->t->top++; /* Push u to TS */
  assert(tmp_G->t->top < 1+tmp_G->v+tmp_G->e && "tmp_G->t->top overflows!");

  if(tmp_G->t->d[u] == tmp_G->t->low[u]){
    (*tmp_G->t->sccC)++; /* Count another SCC */

    /* Rewrite low as identifier */
    *tmp_G->t->SCCsize = 0;
    tmp_G->t->top--;
    assert(0 <= tmp_G->t->top && "tmp_G->t->top below 0");
    do {
      if(1 == tmp_G->t->TSbool[tmp_G->t->TS[tmp_G->t->top]]){
        /* Found a finishing item */
        tmp_G->t->low[tmp_G->t->TS[tmp_G->t->top]] = 1+u;  /* Avoid 0x */
        (*tmp_G->t->SCCsize)++;
        *tmp_G->t->orderP = tmp_G->t->TS[tmp_G->t->top];
        tmp_G->t->orderP++;
        tmp_G->t->TSbool[tmp_G->t->TS[tmp_G->t->top]] = 0; /* Really pop */
      }
      tmp_G->t->top--;
      assert(0 <= tmp_G->t->top && "tmp_G->t->top below 0");
    } while(tmp_G->t->TS[tmp_G->t->top] != u);

    tmp_G->t->SCCsize++; /* Move to the next SCC */
  }
  *tmp_G->t->SCCsize = -1; /* for debug loop */
}

graph
sccRestrict(graph Gin,
	    int *order,
	    int *SCCSA,
      int *SCCED,
      int inPlace
	    )
{
  graph G = Gin;
  // fprintf(stderr, "Before tarjanAllocG\n");
  tarjanAllocG(G, inPlace); // first time: alloc, also checks if Graph changes
  // fprintf(stderr, "After tarjanAllocG\n");
  G->t->top = 0;
  G->t->visited = 0;
  if (!inPlace)
  {
    G->t->orderP = &order[1];
    G->t->order = &order[0];
    G->t->sccC = &order[0];
    G->t->SCCsize = SCCSA;
    G->t->SCCsizeS = SCCSA;
    G->t->SCCnbEdges = SCCED;
  }
  else
  {
    SCCSA = G->t->SCCsize = G->t->SCCsizeS;
    SCCED = G->t->SCCnbEdges;
    order = &G->t->order[0];
    G->t->orderP = &order[1];
  }

  order[0] = 0;

  int *p;
  graph rG = G->t->rG;
  if (inPlace)
  {
    tarjanReset(G);
    if (G->t->t != G->start_E)
      G->t->t = G->start_E;
    rG->t = G->t;
    G = rG;
  }

  // fprintf(stderr, "Before memcpy outDeg/inDeg\n");
  if (!inPlace)
  {
    memcpy(rG->deg[out], G->deg[out], G->v*sizeof(int));
    memcpy(rG->deg[in], G->deg[in], G->v*sizeof(int));
  }
  // fprintf(stderr, "After memcpy outDeg/inDeg\n");

  tmp_G = G;
  
  for(int i = 0; i < G->v; i++){
    if(0 == G->t->d[i])
      TarjanVisit(i);
  }
  // p = NULL;
  // *p = 0; // forces the error
  // act.sa_handler = SIG_DFL;
  // sigaction(SIGSEGV, &act, NULL);
  // fprintf(stderr, "After TarjanVisit\n");

  /* Restricted graph */
  // fprintf(stderr, "Before setup adj list\n");
  for(int u = 0; u < G->v; u++){
    p = G->E[u][out];
    while(-1 != *p){
      if(G->t->low[u] != G->t->low[*p]) {
        rG->e--; /* Remove SCC crossing edges */
        rG->deg[out][u]--;
        rG->deg[in][*p]--;
        assert(rG->deg[in][*p] >= 0);
        assert(rG->deg[out][u] >= 0);
      }
      p++;
    }
  }
  // fprintf(stderr, "After setup adj list\n");

  /* Load info into the new graph */
  // int *t = malloc(2*(rG->v+rG->e)*sizeof(int));
  long unsigned int sizeE = 2*(G->v+G->e)*sizeof(int);
  int *oldT = G->t->t;
  if (0 == G->t->tSize || inPlace) {
    G->t->tSize = sizeE;
    G->t->t = malloc(G->t->tSize);
  } else if (G->t->tSize < sizeE) {
    G->t->tSize = sizeE;
    // G->t->t = realloc(G->t->t, G->t->tSize);
  }
  // fprintf(stderr, "Before setup rG adj list\n");
  int *t = G->t->t;
  for(int u = 0; u < G->v; u++) {
    for(enum adjacencies dir = out; dir <= in; dir++) {
      p = G->E[u][dir];
      rG->E[u][dir] = t;
      while(-1 != *p) {
        if(G->t->low[u] == G->t->low[*p]) {
          *t = *p;
          t++;
        }
        p++;
      }
      *t = -1;
      t++;
    }
  }
  if (inPlace)
  {
    if (rG->start_E && rG->start_E != oldT)
      free(rG->start_E);
    free(oldT);
  }
  rG->start_E = G->t->t;
  // fprintf(stderr, "After setup rG adj list\n");

  // fprintf(stderr, "Before countEdgesPerSCC\n");
  countEdgesPerSCC(G, 1); // TODO: used to check if SCC is complete
  // fprintf(stderr, "After countEdgesPerSCC\n");

#ifndef NDEBUG
  // memset(G->t->d, 0, G->v * sizeof(int)); // allow to run again after G is modified
  for (int i = 0; i < rG->v; i++) {
    assert(rG->deg[in][i] == degree(rG, i, in));
    assert(rG->deg[out][i] == degree(rG, i, out));
  }
#endif

  return rG;
}

int
updateRemFromAdjacencyList(graph G, int v, int *visited)
{
  int *p = G->E[v][out];
  int *q = G->E[v][in];
  int *t = visited;

  while (*p != -1) {
    removeOneVertexFromList(v, G->E[*p][in]);
    // printf("remove %d from in of %d\n", v, *p);
    G->deg[in][*p]--;
    assert(0 <= G->deg[in][*p]);
    G->d[*p] = G->policy_func(G->deg[in][*p], G->deg[out][*p]);
    G->e--;
    if (t) {
      *t = *p;
      t++;
    }
    p++;
  }
  while (*q != -1) {
    removeOneVertexFromList(v, G->E[*q][out]);
    // printf("remove %d from out of %d\n", v, *q);
    G->deg[out][*q]--;
    assert(0 <= G->deg[out][*q]);
    G->d[*q] = G->policy_func(G->deg[in][*q], G->deg[out][*q]);
    G->e--;
    if (t) {
      *t = *q;
      t++;
    }
    q++;
  }
  if (t) *t = -1;
  G->deg[out][v] = 0;
  G->deg[in][v] = 0;
  G->E[v][in][0] = -1;
  G->E[v][out][0] = -1;
  G->d[v] = G->policy_func(0, 0);
  return t - visited;
}

// // merges v into u, puts u in fvs and removes it if self-loop detected
// int // THIS OPERATION IS EXPENSIVE!!!
// updateMerge2Verts(graph G, int v, int u, int *fvs)
// {
//   /* TODO */
//   return 0;
// }

void
removeOneVertexFromList(int toRem, int *list)
{
  int *l = list;
  int *p;
  while (*l != -1 && *l != toRem) l++;
  // we reached either toRem or -1
  p = l;
  assert(*l == toRem || *l == -1);
  while (*l != -1) l++;
  if (l-1 != p) *p = *(l-1);
  *(l-1) = -1;
}

int *
trimG(graph G, int *vertList, int *trimmed, int *procVerts)
{
  int *l = vertList;
  int *t = trimmed;
  int *visited = malloc((G->v+1)*sizeof(int));
  *visited = -1;

  while (*l != -1) {
    int v = *l;
    int is_trimmed = 0;
    // IN0;OUT0 rule
    if (0 == G->deg[in][v] || 0 == G->deg[out][v]) {
      *t = v;
      is_trimmed = 1;
      t++;
      assert(t-trimmed < G->v+1);
      // printf("Trimming %d\n", *l);
      int sizeVis = updateRemFromAdjacencyList(G, v, visited);
      assert(sizeVis < G->v);
      if (sizeVis > 0) {
        t = trimG(G, visited, t, procVerts);
      }
      if (procVerts) (*procVerts)++;
    }
    // IN1;OUT1 rule --> may create self-loops (solved elsewhere)
    // TODO: seems expensive
    if (is_trimmed) {
      removeOneVertexFromList(*l, vertList);
    } else l++;
    assert(l-vertList < G->v+1); // allow for the ending -1
  }
  free(visited);
  *t = -1;
#ifndef NDEBUG
  for (int i = 0; i < G->v; i++) {
    assert(G->deg[in][i] == degree(G, i, in));
    assert(G->deg[out][i] == degree(G, i, out));
    assert(G->d[i] == G->policy_func(G->deg[in][i], G->deg[out][i])); // DOES NOT WORK FOR RANDOM POLICY!
  }
#endif
  return t;
}

int
degree(graph G,
       int v, /* The node */
       enum adjacencies dir
       )
{
  int r = 0;

  int *pv = G->E[v][dir];
  while(-1 != *pv){
    r++;
    pv++;
  }

  return r;
}

int *
getBestSolutionBuffer(graph G)
{
  return G->t->L;
}

static int
pcmp_scores(
  const void *a,
  const void *b,
  void *args
) {
  int *scores = (int *)args;
  int v = *(int *)a;
  int u = *(int *)b;
  int sv = scores[v];
  int su = scores[u];
  return sv - su;
}

static int
pcmp_largeLast(
  const void *a,
  const void *b,
  void *args
) {
  graph G = (graph)args;
  int v = *(int *)a;
  int u = *(int *)b;
  int sv = G->deg[in][v] * G->deg[out][v];
  int su = G->deg[in][u] * G->deg[out][u];
  return sv - su;
}

static int
pcmp_largeFirst(
  const void *a,
  const void *b,
  void *args
) {
  graph G = (graph)args;
  int v = *(int *)a;
  int u = *(int *)b;
  int sv = G->deg[in][v] * G->deg[out][v];
  int su = G->deg[in][u] * G->deg[out][u];
  return su - sv;
}

void
sortVertexArrayOnDegProd_largeFirst(graph G,
  int *array,
  int array_size
      )
{
  qsort_r(array, array_size, sizeof(int), pcmp_largeFirst, G);
}

void
sortVertexArrayOnDegProd_largeLast(graph G,
  int *array,
  int array_size
      )
{
  qsort_r(array, array_size, sizeof(int), pcmp_largeLast, G);
}

void
sortArray_scores(int *array,
  int *scores,
  int array_size
      )
{
  qsort_r(array, array_size, sizeof(int), pcmp_scores, scores);
}

#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
