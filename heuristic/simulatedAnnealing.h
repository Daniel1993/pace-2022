#ifndef _SIMULATED_ANNEALING_H
#define _SIMULATED_ANNEALING_H

#include "sa_state.h"
#include "sa_lib.h"

/* Used for defining calibration values */
double
getTemperature(double p, /* probability */
	       int dE /* Energy delta */
	       );


/* Return pointer to the end of the buffer */
/* The array of vertexes A, does not need */
/* to be a DAG. A discarding heuristic is used. */
int *
executeSA(
	sa_state s, /* State struct to use */
	double hot, /* Hot temperature */
	double decay, /* Cold temperature */
	int *L, /* Buffer for storing solution */
	int improveOnlyAfter /* if set to 1, only updates the final solution on exit */
	  );

int *
guessSA(sa_state s,
  int *L,
  int attempts
      );

void
shuffleSA(
    sa_state s, /* State struct to use */
		int *SCCnbEdges,
    int *SCCsize,
    int curSCC,
    int nbVerts
	  );

#endif /* _SIMULATED_ANNEALING_H */
