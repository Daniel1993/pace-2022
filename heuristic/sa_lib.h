#ifndef SA_LIB_H_GUARD
#define SA_LIB_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "sa_state.h"
#include "graph.h"
#include <stdlib.h>
#include <time.h>

#ifndef INIT_OPS // initial ops to be printed in the file
#define INIT_OPS 5000
#endif

#ifndef THRESHOLD_OPS_NOT_IMPROVING // TODO: fine-tune this parameter

#define THRESHOLD_OPS_NOT_IMPROVING 0.001f
#define THRESHOLD_RESET_TEMPERATURE 0.02f
#define MAX_THRESHOLD_RESET_TEMPERATURE 0.2f

// #define THRESHOLD_OPS_NOT_IMPROVING 1L
// #define THRESHOLD_RESET_TEMPERATURE 16L
// #define MAX_THRESHOLD_RESET_TEMPERATURE 2048L

#endif /* OPS_NOT_IMPROVE*/

#define MAX_REORGANIZE_FACTOR   128.0f
#define MIN_REORGANIZE_FACTOR   0.1f
#define REORGANIZE_FACTOR_DECAY 0.4f

extern double thresholdResetTemperature;
extern double thresholdOpsNotImproving;

// typedef struct graph_ *graph;
// typedef struct state_ *sa_state;

typedef struct SA_ *SA_s;
struct SA_ { // TODO: hide the struct
  long int opsNotImproving;
  long int repOps;// = 0; /* When to report information */
  int *maxE; /* Maximum obtained thus far */
  int *localMaxE; /* Maximum obtained thus far (local) */
  int *currMaxE; /* debug */
  double hot;// = 0.10;
  int hotD;// = -5;
  double cold;// = 0.10;
  int coldD;// = -1;
  double decay;
  unsigned char *inMaxDag;

  // break the execution in batches
  int numberOfBatch;
  int *complexSCCs;

  graph G;
  graph rG;
  sa_state *s; // per SCC
  int maxDagSize;
  int accScore;
  int **vertsPerSCC;
  long *budgetPerSCC;
  int currSCC;
  int currSCCsize;
  int *order;
  int *orderTransSCC;
  int *SCCsize;
  int *SCCnbEdges;
  graph *SCCgraphs;
};

typedef struct SA_parameters_ {
  double hot;// = 0.10;
  int hotD;// = -5;
  double cold;// = 0.10;
  int coldD;// = -1;
  double decay;
} SA_parameters_s;

int
SA_init_G(SA_parameters_s, graph);

int
SA_init_F(SA_parameters_s, const char *filename);

int
SA_reset();

int
SA_get_nbVerts();

int
SA_get_nbEdges();

void
SA_run();

int*  // loop with while(-1 != *bestSolution)
SA_getBestSolution();

void
SA_destroy();

void
SA_set_prof_file(FILE *fp);

// void
// SA_printHeader();

// void
// SA_printLine();

struct timespec
SA_getStartTime();

void
SA_printFVS();

#ifdef __cplusplus
}
#endif

#endif /* SA_LIB_H_GUARD */