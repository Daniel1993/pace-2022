#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "splayTree.h"
#include "bitmap.h"

struct node{
  unsigned int s; /* Number of nodes in sub-tree */
  node c[2]; /* child array */
  node *hook; /* How is this node connected? */
};

struct sTree{
  int n; /* sizeof N */
  node N; /* Node array */
  node root; /* Root of the Splay tree */
};

sTree
allocSTree(int n
           )
{
  sTree t = (sTree)calloc(1, sizeof(struct sTree));
  t->N = (node)calloc(n, sizeof(struct node));
  t->n = n;

  return t;
}

void
freeSTree(sTree t
          )
{
  free(t->N);
  free(t);
}

static void
clearTreeR(node n
	   )
{
  if(NULL != n){
    if(n->hook) *(n->hook) = NULL;
    n->hook = NULL;
    clearTreeR(n->c[left]);
    clearTreeR(n->c[right]);
  }
}

/* Turns into an empty tree */
void
clearTree(sTree t
          )
{
  clearTreeR(t->root);
}

static void
clearNode(sTree t,
	  int i
          )
{
  bzero(&t->N[i], sizeof(struct node));
}

int
size(sTree t
     )
{
  int r = 0;
  if(NULL != t->root)
    r = t->root->s;

  return r;
}

// TODO: why the below does not work?
// #define father(_n) (&(_n)[(struct node*)(_n)->hook - (struct node*)(_n)])
static node
father(node n
       )
{
  assert(NULL != n->hook && "Father call on root node.");
  long int d = (long int) n->hook - (long int) n;
  if (0 > d)
  {
    d -= sizeof(struct node);
    d++;
  }
  node f = &n[d/(int)sizeof(struct node)];
  return f;
}

#define childIdx(_n, _f) (_f->c[right] == _n ? right : left)
/* static enum childI
childIdx(node n, // The child
	 node f // The father
	 )
{
  enum childI i = left;

  if(f->c[right] == n)
    i = right;

  return i;
} */

// assumes left==0 and right==1, use the XOR operator to flip the direction
#define other(_dir) (_dir ^ 1)
/* static enum childI
other(enum childI dir
      )
{
  enum childI r = left;
  if(r == dir)
    r = right;

  return r;
} */

#define KEY(_t, _n) ((node)(_n) - (node)(_t)->N)
int
key(sTree t,
    node n
    )
{
  return KEY(t, n);
}

static sTree auxT;
static int *auxL;
static int *auxMap;
static enum childI auxD;

/* // TODO: init this elsewhere (iterative version seems a bit more efficient)
#define MAX_POSSIBLE_VERTS 1024*1024
static node tmp_node_stack[MAX_POSSIBLE_VERTS];
static node *tmp_node_stack_f;
static void
getInorderR_mapNULL_it(node n
            )
{
  node tmpN = n;
  node auxN;
  *(tmp_node_stack_f++) = tmpN;
  while (NULL != (auxN = tmpN->c[auxD]))
  {
    *(tmp_node_stack_f++) = auxN;
    tmpN = auxN;
  }

  assert(tmp_node_stack_f < tmp_node_stack+MAX_POSSIBLE_VERTS && "Not enough space in the stack!");

  while (tmp_node_stack_f-- > tmp_node_stack)
  {
    tmpN = *tmp_node_stack_f;
    *(auxL++) = KEY(auxT, tmpN);
    if(NULL != (auxN = tmpN->c[other(auxD)]))
    {
      *(tmp_node_stack_f++) = auxN;
      tmpN = auxN;
      while (NULL != (auxN = tmpN->c[auxD]))
      {
        *(tmp_node_stack_f++) = auxN;
        tmpN = auxN;
      }
    }
  }
} */

static void
getInorderR_mapNULL(node n
            )
{
  if(NULL != n->c[auxD])
    getInorderR_mapNULL(n->c[auxD]);

  *(auxL++) = KEY(auxT, n);

  if(NULL != n->c[other(auxD)])
    getInorderR_mapNULL(n->c[other(auxD)]);
}

static void
getInorderR(node n
            )
{
  if(NULL != n->c[auxD])
    getInorderR(n->c[auxD]);

  *(auxL++) = auxMap[KEY(auxT, n)];

  if(NULL != n->c[other(auxD)])
    getInorderR(n->c[other(auxD)]);
}

/* static int *
getInorderR(sTree t,
            node n,
            int *L,
            int *map,
            enum childI d
            )
{
  if(NULL != n->c[d])
    L = getInorderR(t, n->c[d], L, map, d);

  if(NULL == map)
    *L = KEY(t, n);
  else
    *L = map[KEY(t, n)];
  L++;

  if(NULL != n->c[other(d)])
    L = getInorderR(t, n->c[other(d)], L, map, d);

  return L;
} */

int *
getInorder(sTree t,
           int *L,
	         int *map,
	         enum childI d
          )
{
  int *F = L;

  if(NULL != t->root)
  {
    // F = getInorderR(t, t->root, L, map, d);
    auxT = t;
    auxL = L;
    auxMap = map;
    auxD = d;
    // tmp_node_stack_f = &(tmp_node_stack[0]);
    if (NULL == map)
      getInorderR_mapNULL(t->root); //getInorderR_mapNULL_it(t->root);
    else
      getInorderR(t->root);
    F = auxL;
  }

  return F;
}

#define resetNode(_n) \
  *((_n)->hook) = (_n); \
  (_n)->s = 1; \
  if(NULL != (_n)->c[left]) \
    (_n)->s += (_n)->c[left]->s; \
  if(NULL != (_n)->c[right]) \
    (_n)->s += (_n)->c[right]->s; \
//
/* static void
resetNode(node n
          )
{
  *(n->hook) = n;

  n->s = 1;
  if(NULL != n->c[left])
    n->s += n->c[left]->s;
  if(NULL != n->c[right])
    n->s += n->c[right]->s;
} */

#ifndef NDEBUG
static void
checkTree(sTree t
          )
{
#if 0
  for(int i = 0; i < t->n/*(int)(sizeof(t->N)/sizeof(struct node))*/; i++){
    if(NULL != t->N[i].hook) {
      // printf("check %p father (%p)\n", &(t->N[i]), father(&(t->N[i])));
      assert(&(t->N[i]) != father(&(t->N[i])) && "Found loop in tree");
    }
  }
#endif
}
#endif /* NDEBUG */

static void
splay(sTree t,
      node n
      )
{ /* WARNING: Apply only to nodes in the tree */
  assert(NULL != n->hook && "Splay call on out of tree node");
#ifndef NDEBUG
  checkTree(t);
#endif /* NDEBUG */

  while (&(t->root) != n->hook)
  {
    node f = father(n); /* Father node */
    assert(f != n && "Strange loop in tree");
    enum childI fi = childIdx(n, f);
    if (&(t->root) == f->hook)
    {
      /* Simple Zig transformation
         |
         |     f    =>    n
         |    / \   =>   / \
         |   n   C  =>  A   f
         |  / \     =>     / \
         | A   B    =>    B   C
         |
         | fi = /   ;  other(fi) = \   */

      n->hook = f->hook; /* Copy hook */
      f->c[fi] = n->c[other(fi)];
      if(NULL != f->c[fi])
        f->c[fi]->hook = &f->c[fi];
      f->hook = &n->c[other(fi)];
    } else {
      node gf = father(f); /* Grand Father node */
      assert(gf != n && gf != f && "Weird loop in tree");
      enum childI gfi = childIdx(f, gf);

      if (fi == gfi) /* Zig Zig */
      {
        /* Zig Zig transformation
          |
          |       gf   =>    n
          |      /  \  =>   / \
          |     f    D =>  A   f
          |    / \     =>     / \
          |   n   C    =>    B   gf
          |  / \       =>       /  \
          | A   B      =>      C    D
          |
          | fi = /   ;  other(fi) = \    */

        n->hook = gf->hook;
        f->hook = &n->c[other(fi)];
        gf->hook = &f->c[other(fi)];

        f->c[fi] = n->c[other(fi)];
        if (NULL != f->c[fi])
          f->c[fi]->hook = &f->c[fi];
        gf->c[fi] = f->c[other(fi)];
        if (NULL != gf->c[fi])
          gf->c[fi]->hook = &gf->c[fi];
      }
      else /* Zig Zag */
      {
        /* Zig Zag transformation
          |
          |     gf    =>       n
          |    /  \   =>      / \
          |   f    D  =>     /   \
          |  / \      =>    /     \
          | A   n     =>   f       gf
          |    / \    =>  / \     /  \
          |   B   C   => A   B   C    D
          |
          | fi = \   ;  gfi = /     */

        n->hook = gf->hook;
        f->hook = &n->c[gfi];
        gf->hook = &n->c[fi];

        f->c[fi] = n->c[gfi];
        if(NULL != f->c[fi])
          f->c[fi]->hook = &f->c[fi];
        gf->c[gfi] = n->c[fi];
        if(NULL != gf->c[gfi])
          gf->c[gfi]->hook = &gf->c[gfi];
      }
      resetNode(gf);
    }
    resetNode(f);
    resetNode(n);
  }
  resetNode(n);

#ifndef NDEBUG
  checkTree(t);
#endif /* NDEBUG */
}

/* Return value stored in this node. */
int
value(sTree t,
      node n,
      int *cv /* complementar value */
      )
{ /* WARNING: Apply only to nodes in the tree */
  assert(NULL != n->hook && "Splay call on out of tree node");

  int r = 0;
  splay(t, n);
  if(NULL != n->c[left])
    r = n->c[left]->s;

  if(NULL != cv){
    *cv = 0;
    if(NULL != n->c[right])
      *cv = n->c[right]->s;
  }

  return r;
}

/* Return node pointer by index */
node
getNode(sTree t,
        int i
        )
{
  node r = NULL;
  if(NULL != t->N[i].hook)
    r = &(t->N[i]);

  return r;
}

node
getNodeInOrder(sTree t,
	       int i
	       )
{
  node r = t->root;

  assert(NULL != r && "Request node on empty tree.");

  int l = 0; /* Number of elements on the left */

  while(NULL != r && l != i){
    if(NULL != r->c[left])
      l += r->c[left]->s;

    if(l != i){
      enum childI d = left;
      if(l < i)
	      d = right;

      r = r->c[d];
    }
  }

  return r;
}

void
insertKey(sTree t,
	  int k
	  )
{
  assert(NULL == t->N[k].hook && "Inserting existing node");
  // if (NULL != t->N[k].hook) { printf("duplicate!\n"); return; }


  if(NULL == t->root){ /* Insert on empty tree */
    reRoot(t, k, right);
  } else {
    node floor;
    node ceil;
    roundSt(t, k, &floor, &ceil);

    if(NULL == floor)
      reRoot(t, k, right);
    else if(NULL == ceil)
      reRoot(t, k, left);
    else
      insertN(t, floor, k, right);
  }
}

void
insertInorderKey(sTree t,
		 int p, /* position */
		 int k
		 )
{
  assert(NULL == t->N[k].hook && "Inserting existing node");

  node *h = &t->root;
  node u = *h;

  while(NULL != u){
    int lsize = 0;
    if(NULL != u->c[left])
      lsize = u->c[left]->s;

    if(lsize < p){
      p -= 1 + lsize;
      h = &u->c[right];
    } else {
      h = &u->c[left];
    }
    u = *h;
  }

  clearNode(t, k);
  node v = &t->N[k];
  v->hook = h;
  *h = v;
  splay(t, v);
}


/* Insert this node into the tree */
void
insertN(sTree t,
        node n, /* Reference node */
        int i, /* Node to insert */
        enum childI d /* Direction of reference */
        )
{ /* WARNING: Apply only to nodes in the tree */
  assert(NULL != n->hook && "Splay call on out of tree node");
  assert(NULL == t->N[i].hook && "Inserting existing node");

  if(NULL != n->c[d]){
    n = n->c[d];
    d = other(d);
    while(NULL != n->c[d])
      n = n->c[d];
  }

  clearNode(t, i);
  node v = &t->N[i];
  v->hook = &n->c[d];
  *(v->hook) = v;
  splay(t, v);
}

/* Make the node corresponding to i the new root */
void
reRoot(sTree t,
       int i, /* New root index */
       enum childI d /* Direction of reference */
       )
{
  assert(NULL == t->N[i].hook && "Rooting to active node");

  clearNode(t, i);

  t->N[i].c[d] = t->root;
  if(NULL != t->N[i].c[d])
    t->N[i].c[d]->hook = &(t->N[i].c[d]);

  t->N[i].hook = &t->root;
  resetNode(&t->N[i]);
}

void
removeN(sTree t,
        node n /* Node to remove */
        )
{ /* WARNING: Apply only to nodes in the tree */
  assert(NULL != n->hook && "Splay call on out of tree node");

  if(1 == t->root->s) /* Tree has only 1 node */
    t->root = NULL;
  else {  /* Tree has several nodes */
    splay(t, n);

    node o = n->c[left];
    if(NULL == o)
      o = n->c[right];
    else if(NULL != n->c[right]) { /* General case */
      while(NULL != o->c[right])
        o = o->c[right];

      splay(t, o); /* Splay the other node */
      splay(t, n); /* Splay n again */

      /* Removal transform
        |
        |     n    =>
        |    / \   =>
        |   o   B  =>    o
        |  /       =>   / \
        | A        =>  A   B
        |
        | d = \     ; other(d) = /     */

      o->c[right] = n->c[right];
      o->c[right]->hook = &o->c[right];
    }

    o->hook = &t->root;
    resetNode(o);
  }
  n->hook = NULL;
}

void
splitSt(sTree t,
	int k, /* The key */
	enum childI d /* Direction */
	)
{
  node floor;
  node ceil;

  roundSt(t, k, &floor, &ceil);

  switch(d){
  case left:
    if(NULL == floor)
      clearTree(t);
    else{
      splay(t, floor);
      clearTreeR(floor->c[right]);
    }
    break;

  case right:
    if(NULL == ceil)
      clearTree(t);
    else{
      splay(t, ceil);
      clearTreeR(ceil->c[left]);
    }
    break;

  default:
    assert(0 && "Invalid split direction ");
    break;
  }
}

/* Extreme movement */
static node
extreme(sTree t,
	node n,
	enum childI d
	)
{
  assert(NULL != n && "extreme call on NULL node.");

  while(NULL != n->c[d])
    n = n->c[d];

  splay(t, n);

  return n;
}

node
diressor(sTree t,
	 node n, /* Reference node */
	 enum childI d
	 )
{
  assert(NULL != n && "diressor call on NULL node.");

  splay(t, n); /* Make it into root */
  if(NULL != n->c[d]){
    n = n->c[d];

    n = extreme(t, n, other(d));
  }

  return n;
}

node
dirum(sTree t,
      enum childI d
      )
{
  node r = t->root;

  if(NULL != r)
    r = extreme(t, r, d);

  return r;
}

void
roundSt(sTree t,
        int k, /* The key */
	node *floor,
	node *ceil
        )
{
  node n = t->root;
  *floor = NULL;
  *ceil = NULL;
  node last = NULL;

  while(NULL != n && KEY(t, n) != k){
    last = n;

    if(KEY(t, n) < k){
      *floor = n;
      n = n->c[right];
    } else {
      *ceil = n;
      n = n->c[left];
    }
  }

  if(KEY(t, n) == k){
    last = n;
    *floor = n;
    *ceil = n;
  }

  if(NULL != last)
    splay(t, last);
}
