#include <assert.h>
#include <string.h>
#include <stdlib.h>

// #include <bsd/stdlib.h>
#include "WELL512a.h"
#include "queue.h"
#include "splayTree.h"
#include "reorganize.h"
#include "sa_state.h"

struct organizer{
  graph G;
  char *A; /* Break node booleans */
  int *R; /* Restart order */

  unsigned int *d; /* discovery time */
  unsigned int *low; /* low for Tarjan like analysis */

  int *L; /* Tarjan Stack */
  /* char *Lb; /\* Stack booleans *\/ */
  int *d2L; /* Discovery to stack position */
  int top;

  sTree t; /* Resulting order */
  unsigned cnt; /* Heuristic counter */
};

unsigned int time; /* Discovery times */

static int
hasDQ(organizer o, unsigned int d)
{
  int r = (unsigned int) -1 != d;
  r = r && o->d2L[d] <= o->top;
  r = r && o->d[o->L[o->d2L[d]]] == d;

  return r;
}

static int
inStack(organizer o, int u)
{
  int r = (unsigned int) -1 != o->d[u];
  r = r && o->d2L[o->d[u]] <= o->top;
  r = r && u == o->L[o->d2L[o->d[u]]];

  return r;
}

static void
popOS(organizer o)
{
  o->top--;

  /* fprintf(stderr, "/ "); */
}

static void
pushOS(organizer o, int u)
{
  o->top++;

  /* fprintf(stderr, "* "); */

  o->d2L[o->d[u]] = o->top;
  o->L[o->top] = u;
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

/* Creates a new organizing structure */
organizer
allocO(graph G,
       sTree t
       )
{
  organizer o = malloc(sizeof(struct organizer));

  o->G = G;
  o->A = calloc(G->v, sizeof(char));
  o->R = malloc(G->v*sizeof(int));
  for(int i = 0; i < G->v; i++)
    o->R[i] = i;

  o->d = malloc(G->v*sizeof(int));
  o->low = malloc(G->v*sizeof(int));

  o->top = -1;
  o->L = malloc((G->v+G->e)*sizeof(int));
  o->d2L = malloc(G->v*sizeof(int));

  o->t = t;
  o->cnt = 0;

  return o;
}

void
freeO(organizer o
      )
{
  free(o->A);
  free(o->R);

  free(o->d);
  free(o->low);

  free(o->L);
  free(o->d2L);

  o->t = NULL;

  free(o);
}

// TODO: place it in a better place
extern float reorganize_factor;

int
checkO(organizer o,
       int dE
       )
{
  // TODO: set this factor to balance calls to dfs
  // with algorithm precision
  int r = o->cnt > reorganize_factor * (o->G->v + o->G->e);

  if(!r){
    o->cnt += dE;
  }

  return r;
}

/* Using global vars to optimize stack space */
static organizer Go;
static node tu; /* Insert after this node */
static int Br; /* Current break node under eval */

/* Shuffle and put break nodes at the end */
static void
breakPartition(int n,
	       int* Adj,
	       char *A /* Active booleans */
	       )
{
  /* Shuffle */
  for(int i = 0; i < n-1; i++){
    int rd = WELLRNG512a_u4()%(n-i);
    swapI(&Adj[rd], &Adj[i]);
  }
  /* Partition, put break nodes at the end */
  int i = 0;
  while(i < n){
    while(i < n && 0 == A[Adj[i]]) i++;
    while(i < n && 0 != A[Adj[n-1]]) n--;
    if(i < n){
      n--;
      swapI(&Adj[n], &Adj[i]);
      i++;
    }
  }
}

/* Double call solves all your problems */
static void
dfs(int v /* Vertex to start from */
    )
{
  organizer o = Go;

  /* fprintf(stderr, "dfs_visit_start(%d)\n", v); */

  if(-1 == Br && 1 == o->A[v]){ /* add test node */
    Br = v;
    o->A[v] = 0;
    /* fprintf(stderr, "T"); */
  }

  if(0 == o->A[v]){ /* Visitable */
    time++;
    o->d[v] = time; /* Do not init low to exclude self-loops */
    pushOS(o, v);

    /* Recurse */
    int *Adj = o->G->E[v][out];
    int n = o->G->deg[out][v];
    // int *tj = Adj;

    // int n = 0;
    // while(-1 != *Adj){ n++; Adj++; }

    // Adj = tj;
    breakPartition(n, Adj, o->A);

    while (-1 != *Adj) {
      int u = *Adj;
      Adj++;
      if((unsigned int)-1 == o->d[u])
        dfs(u);

      if(-1 != Br){ /* Only when you are @ or bellow test node */
        /* Is u in the stack? */
        if(o->low[v] > o->d[u] &&
            inStack(o, u))
          o->low[v] = o->d[u];

        /* Is the node with d=low[u] in stack? */
        if(o->low[v] > o->low[u] &&
          hasDQ(o, o->low[u]))
          o->low[v] = o->low[u];
      }
    }

    if(v == Br){ /* Judge node */
      if(o->low[v] <= o->d[v]){ /* Cycle */
        o->A[v] = 2;
        do
          popOS(o);
        while(v != o->L[o->top+1]);
      }
      Br = -1; /* Ready for new tests */
    }

    if(o->low[v] > o->d[v]){ /* No cycle so just pop */
      assert(v == o->L[o->top] && "Single SCC invariant.");
      popOS(o);  /* Otherwise wait */
    }

    if(NULL != tu && 0 == o->A[v]){
      /* fprintf(stderr, "A:%d ", v); */
      insertN(o->t, tu, v, right);
    }
  } else
    o->A[v] = 2; /* Make it non-insertable */

  /* fprintf(stderr, "d:%d l:%d %d", o->d[v], o->low[v], v); */
  /* if(2 == o->A[v]) */
  /*   fprintf(stderr, "X"); */

  /* fprintf(stderr, ") "); */

  /* fprintf(stderr, "dfs_visit_end(%d)\n", v); */
}

/* Inserts topological order into t */
int flag_did_reorganize = 0;
int flag_reorganize_improve = 0;
void
reorganize(organizer o,
	   int *ln, /* list size */
	   int *idx, /* idx in the list */
	   int *L /* List of break nodes */
	   )
{
  o->cnt = 0;
  flag_did_reorganize = 1;

  for(int i = 0; -1 != L[i]; i++)
    o->A[L[i]] = 1;

  clearTree(o->t);
  reRoot(o->t, o->G->v, left);

  /* First DFS establishes nodes */
  /* for(int i = 0; i < o->G->v && 0 == o->A[Rst[i]]; i++) */
  /*   fprintf(stderr, "B:%d ", Rst[i]); */
  /* fprintf(stderr, "\n"); */

  Go = o;

  Br = -1; /* Ready for testing */
  tu = NULL;
  time = -1;
  for(int i = 0; i < o->G->v; i++){
    o->d[i] = time;
    o->low[i] = time;
  }

  breakPartition(o->G->v, o->R, o->A);
  for(int i = 0; i < o->G->v; i++)
    if((unsigned int)-1 == o->d[o->R[i]]){
      dfs(o->R[i]);
      assert(-1 == o->top && "Messing up the stack.");
    }

  /* Second DFS gets order */
  Br = o->G->v+1; /* No tests */
  tu = getNode(o->t, o->G->v); /* Just collect */
  time = -1;
  for(int i = 0; i < o->G->v; i++){
    o->d[i] = time;
    o->low[i] = time;
  }

  breakPartition(o->G->v, o->R, o->A);
  for(int i = 0; i < o->G->v && 0 == o->A[o->R[i]]; i++)
    if((unsigned int)-1 == o->d[o->R[i]]){
      dfs(o->R[i]);
      assert(-1 == o->top && "Messing up the stack.");
    }

  for(int i = 0; -1 != L[i]; i++){
    if(0 == o->A[L[i]]){ /* Remove from candidate */
      (*ln)--;
      idx[L[*ln]] = i;
      idx[L[i]] = -1;
      L[i] = L[*ln];
      L[*ln] = -1;
      i--;
      flag_reorganize_improve = 1;
    } else
      o->A[L[i]] = 0;
  }

  /* Remove sentinel from the tree */
  removeN(o->t, getNode(o->t, o->G->v));
}
