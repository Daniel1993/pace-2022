#ifndef _SA_STATE_H
#define _SA_STATE_H

#include "graph.h"

#ifdef POSTSCRIPT
extern FILE *LOG;
extern int POST_NB_VERT;
#endif /* POSTSCRIPT */

typedef struct sa_state_ *sa_state;

sa_state
allocS(graph G
       );

void
freeS(sa_state s
      );

/* Set sa_state vertex to s */
void
prepareS(sa_state s,
	 int* A, /* Array with vertexes */
	 int l	 /* length of the array */
	 );
void
prepareS_vertCover(sa_state s,
	 int* A, /* Array with vertexes */
	 int l	 /* length of the array */
	 );
void
prepareS_greedy(sa_state s,
         int*  A, /* Array with vertexes */
         int   l	/* length of the array */
);
void
prepareS_empty(sa_state s
	 );

int
getNbCandidates(sa_state s
      );

void
tryAllCandidates(sa_state s,
  int pdE
      );

void
shuffleCandidates(sa_state s
      );

int
forceReorganize(sa_state s
      );

int
tryReorganize(sa_state s
      );

void
decayAndReorganize(sa_state s
      );

/* Chooses a potential transition */
void
choose(sa_state s
       );

/* Returns the proposed energy variation */
int
check(sa_state s,
      int dE
      );

/* Executes validated transition */
/* Returns number of new candidates */
void
execute(sa_state s
       );

/* Return sa_state energy */
int
getE(sa_state s
     );

void
printS(sa_state s
       );
void
printS_v2(sa_state s);

/* Load sa_state into list L */
int *
getSketch(sa_state s,
	  int *L
          );

#endif /* _SA_STATE_H */
