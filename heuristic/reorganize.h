#ifndef _REORGANIZE_H
#define _REORGANIZE_H

#include "graph.h"
#include "sa_state.h"

typedef struct organizer *organizer;

/* Creates a new organizing structure */
organizer
allocO(graph G,
       sTree t
       );
       
void
freeO(organizer o
      );

/* Checks to see if the heuristic amortizes.
   Otherwise increment internal counter. */
int
checkO(organizer o,
       int dE
       );

/* Inserts topological order into t */
void
reorganize(organizer o,
	   int *ln,
          int *idx, /* idx in the list */
	   int *L /* List of break nodes */
	   );

#endif /* _REORGANIZE_H */
