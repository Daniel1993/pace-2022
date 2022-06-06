#ifndef _GRAPH_H
#define _GRAPH_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

// NOTE: in C++ struct graph also defines a keyword "graph"
typedef struct graph_ *graph;
typedef struct TarjanMetadata_ *TarjanMetadata;

enum adjacencies {out = 0, in = 1};
#define OTHER_ADJ(_dir) ((_dir)^1)

// TODO: hide these 2 (now they are used somewhere else)
struct graph_ {
  long v; /* Number of vertexes */
  long e; /* Number of edges */
  int *deg[2]; // keep track of the in-/out-degree count for each vertex

  // array to keep the score of the vertexes (not used by all solutions)
  // can also be the translation table to the vertices in some SCC in other graph
  int *d;
  int *score; // TODO

  // TODO
  // int *dU; // ranking having in account unique vertices per cycle

  // outputs the score for that vertex (e.g., product of inDeg and outDeg, not used in all solutions)
  int (*policy_func)(int inDeg, int outDeg);
  int *(*E)[2]; /* Array of pairs of pointers */
  int *(*Eaux)[2]; /* used to keep info per edge */
  int *start_E;
  int *start_Eaux;
  int *inFVS; // due to contractions
  int inFVSs;
  TarjanMetadata t;
};

struct TarjanMetadata_ {
  graph rG;
  char *TSbool;
  int *TS;
  int top;
  int *low;
  int *d;
  int *t;
  int *L;
  size_t tSize;
  int visited;
  int *order; // ref to begining of orderP
  int *orderP;
  int *sccC; /* SCC counter */
  int *SCCsize;
  int *SCCsizeS; // start of SCCsize
  int *SCCnbEdges;
};


void
sortVertexArrayOnDegProd_largeFirst(graph G,
  int *array,
  int array_size
      );

void
sortVertexArrayOnDegProd_largeLast(graph G,
  int *array,
  int array_size
      );

void
sortArray_scores(int *array,
  int *scores,
  int array_size
      );

graph
allocG(long v, long e);

void
scoreVertexByLoopSize(graph G);

graph
extractSCC(graph G, int *sccList, int lenScc);

void
freeG(graph G);

// matrix accessed via idx_row*nbVertices+idx_col
graph
fromSquareMat(long nbVertices, unsigned char *mat);

graph
loadG(FILE *stream);

void
freeG(graph G);

void
printG(graph G);

int
contractG(graph G
	    );

/* Return a new graph without edges that cross SCCs. Also returns the
   component graph in reverse topological order and each SCC in topological
   order (approximately). */
graph
sccRestrict(graph G,
	    int *order,
	    int *SCCsize,
	    int *SCCnbEdges,
      int inPlace
	    );

void
setEaux(graph G);

void // TODO: move to some util file
removeOneVertexFromList(int toRem, int *list);

int
updateRemFromAdjacencyList(graph G, int v, int *visited);

int *
trimG(graph G, int *vertList, int *trimmed, int *procVerts);

/* Return in or out degree */
int
degree(graph G,
       int v, /* The node */
       enum adjacencies dir
       );

int *
getBestSolutionBuffer(graph G);

#ifdef __cplusplus
}
#endif

#endif /* _GRAPH_H */
