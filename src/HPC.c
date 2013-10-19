/* Hyperset Classification
 *
 * Summary: This program, when compiled with the flags
 *              gcc -std=c99 -DLEVEL=x HPC.c
 *          or by using the included Makefile, will attempt to
 *          determine the number of unique minimal hypersets
 *          within the set of all possible APGs of level x
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>

#define U_INT64 unsigned long long
#define U_INT32 unsigned long
#define GRAPH_T unsigned char
#define ONE8 1u
#define ONE32 1ul
#define ONE64 1ull

// define default bounds if not specified
#ifndef CMIN
#define CMIN ONE64<<(LEVEL*LEVEL-2)
#endif
#ifndef CMAX
#define CMAX (ONE64<<(LEVEL*LEVEL))-1
#endif

// function prototypes
bool is_apg(GRAPH_T *graph);
bool is_minimized(GRAPH_T *graph);
void separate_p (GRAPH_T *graph, GRAPH_T *partitions, U_INT32 p, U_INT32 num_partitions);
bool can_separate(GRAPH_T *graph, GRAPH_T a, GRAPH_T b, GRAPH_T *partitions, U_INT32 num_partitions);
bool is_unique(GRAPH_T *graph, U_INT64 orig);
U_INT64 next_permutation(GRAPH_T *graph, GRAPH_T *indexes);
void swap_rowcol(GRAPH_T *graph, GRAPH_T *indexes, GRAPH_T x, GRAPH_T y);
U_INT64 graph2int(GRAPH_T *graph);
static void *doit( void *arg );

// thread-specific values passed from main
typedef struct {
  pthread_t tid;
  U_INT64 cmin;  // min candidate for a given LEVEL graph
  U_INT64 cmax;  // max candidate for a given LEVEL graph
  U_INT64 count;
  bool done;
} bounds;

int main(void)
{
    U_INT64 cmin = CMIN; // min candidate for a given LEVEL graph
    U_INT64 cmax = CMAX; // max candidate for a given LEVEL graph

    size_t nBounds = 2;
    bounds* boundsArray = (bounds*) calloc( nBounds, sizeof(bounds) );

    bounds b1 = boundsArray[0];
    bounds b2 = boundsArray[1];

    b1.count = 0;
    b1.done = false;
    b1.cmin = cmin;
    // because we filter out "chunks" of candidates who don't end in a sink,
    // we need to make sure that b1's max, and b2's min are not in the middle
    // of a "chunk"
    b1.cmax = ((((cmax - cmin) / 2u) + cmin) & ~((ONE64<<LEVEL)-1)) -1;

    b2.count = 0;
    b2.done = false;
    b2.cmin = b1.cmax + 1;
    b2.cmax = cmax;

#ifdef DEBUG
    printf("B1: min %llu \t max %llu\nB2: min %llu \t max %llu\n",
      b1.cmin,b1.cmax,b2.cmin,b2.cmax);
#endif

    pthread_create( &b1.tid, NULL, &doit, &b1 );
    pthread_create( &b2.tid, NULL, &doit, &b2 );

#ifdef SECTION
    printf("Profiling section %i (MIN: %llu, MAX: %llu)of level %i...\n",
      SECTION,CMIN,CMAX,LEVEL);
#else
    printf("Profiling level %i ...\n",LEVEL);
#endif

    while (1)
      if (b1.done && b2.done)
        break;
      else
#if LEVEL==7
        sleep(60);
#else
        sleep(1);
#endif

    U_INT64 total_count = b1.count + b2.count;
    free(boundsArray);
    printf("\nLevel %i has %llu\n",LEVEL,total_count);

    return 0;
}

// threads created by main() that do the actual work start here
static void *doit( void *arg )
{
    bounds *bnd = (bounds *) arg;
    pthread_detach( pthread_self() );

    U_INT64 candidate;      // binary str represents a square adjacency matrix
    U_INT64 mask;           // binary mask
    U_INT32 row;            // iterator
    U_INT64 my_candidate, temp;
    GRAPH_T graph[LEVEL] = {0}; // graph represented as an array
    mask = (ONE64<<LEVEL)-1;    // mask of LEVEL number of ones to isolate rows
#ifdef DEBUG
    printf("MIN: %llu \nMAX: %llu \nmask: %llu\n",cmin,cmax,mask);
#endif

  candidate = bnd->cmin;
  while (1)
  {
    if (candidate > bnd->cmax) break;
    my_candidate=temp=candidate;

#ifdef DEBUG
    printf("\nCandidate: %llu:\n",temp);
#endif

    candidate += mask+1;

    // create graph from candidate
    for (row=0; row<LEVEL; ++row)
    {
      graph[LEVEL-1-row] = (GRAPH_T)(temp & mask);
      temp >>= LEVEL;
    }

    // check for all properties
    if (is_apg(graph) && is_minimized(graph)
       && is_unique(graph, my_candidate))
    {
#ifdef SHOWNUM
        printf("\n%llu",pay.candidate-1);
#endif
#ifdef DEBUG
        printf("\tcounted");
#endif
        ++bnd->count;
    }
  }

     bnd->done = true;

   return NULL;
}

// check to see if all nodes can be traversed.
// also filter out many candidates if sinks are found in places
// other than the last row.
// also filter for candidates that aren't traversed in the order
// their of their rows (throws out many duplicates, but not all)
bool is_apg(GRAPH_T *graph)
{
  U_INT32 row = 0;  // iterator
  U_INT32 mask = 1; // binary mask
  U_INT32 visited = ONE32<<(LEVEL-1); // mark row 0 as visited

#ifdef DEBUG
  printf("\n Check if APG");
#endif
  mask = ONE32<<(LEVEL-1); // start mask at column 0 (edge to root)

  // check all but last row
  for(row=0; row<(LEVEL-1); ++row)
  {
    // By the time each row is iterated to, it should have been
    // visited (as well as all rows before it), otherwise, even
    // if it's an APG, it's not a candidate, because I choose to
    // only allow candidates that are visited in order.

    // if this row has been visited and is not a sink
    if ((visited & mask) && graph[row])
      visited |= graph[row]; // add this row's edges to visited list
    else
    {
#ifdef DEBUG
      printf("\tcan't traverse in order");
#endif
      return false;
    }
    mask >>= 1;
  }

  // check if last row is visited (we ensure it is a sink in the main loop)
  if (visited & mask)
  {
    GRAPH_T ones = 0, zeros = 0;
#ifdef DEBUG
    printf("\tis APG");
#endif
    // now check if zero before a 1, if so, this isn't a valid candidate
    mask = ONE8<<(LEVEL-2);
    while (!(mask & ONE8))
    {
      if (mask & graph[0])
      {
        ++ones;
        if (zeros) return false;
      }
      else
        ++zeros;

      mask >>= 1;
    }
    return true;
  }

#ifdef DEBUG
    printf("\tnot an APG");
#endif
  return false;
}

// attempt to minimize the graph
bool is_minimized(GRAPH_T *graph)
{
#ifdef DEBUG
  printf("\n Check if minimized");
#endif
#if LEVEL>2
  U_INT32 row=0;            // iterator
  U_INT32 num_partitions=2; // at least two partitions(sink and root)
  U_INT32 p;                // index of current partition
  GRAPH_T partitions[LEVEL] = {0}; // room to separate all nodes
  GRAPH_T sinkmask = 1;     // binary mask to check for edge to sink
  GRAPH_T mask;             // binary mask
  bool modified = false;    // for checking if P(k)==P(k-1)

  mask = ONE8<<(LEVEL-1);   // mask positioned under root

  // separate partitions that point to sink and those that don't
  partitions[0]=1; // place sink in first partition
  for (row=0; row<LEVEL-1; ++row)
  {
    partitions[(sinkmask & graph[row])?1:2] |= mask;
    mask >>= 1;
  }
  if (partitions[2]) num_partitions = 3;

  do
  {
    modified=false;
    p=1;
    while (p<num_partitions && num_partitions<LEVEL)
    {
      // attempt to separate partition p
      separate_p(graph,partitions,p,num_partitions);
      if(partitions[num_partitions])
      {
        ++num_partitions;
        modified = true; // partitions have changed
      }
      ++p;
    }
  }
  while (modified); // repeat while P(k) != P(k-1)

  // by now, all partitions are as separate as possible
  if (num_partitions == LEVEL)
    return true; // minimized
  else
    return false;
#endif  // end LEVEL>2
  return true;
}

// helper function for splitting up a single partition
void separate_p (GRAPH_T *graph, GRAPH_T *partitions, U_INT32 p,
                U_INT32 num_partitions)
{
  GRAPH_T a=0, b;           // index of two nodes we are comparing
  GRAPH_T mask;             // binary mask for locating nodes
  GRAPH_T separated=0;      // remember nodes that separate from 'a'
  mask = ONE8<<(LEVEL-1); // position amask to bit for row 0

  // find the first node in this partition
  while(mask>2) // don't bother with sink or if only 1 non-sink
  {
    if(partitions[p] & mask)
      break;
    ++a;
    mask >>= 1;
  }
  if (mask<4) return; // partition empty or only contains sink

  // check next potential
  b=a+1;
  mask >>= 1;

  while (mask>1) // second node won't be the sink
  {
    if(partitions[p] & mask) // found a second node to compare with 'a'
    {
      if (can_separate(graph,a,b,partitions,num_partitions) ||
          can_separate(graph,b,a,partitions,num_partitions))
        separated |= mask;
    }
    ++b;
    mask >>= 1;
  }

  partitions[num_partitions]=separated;// place separated nodes in new partition
  partitions[p] ^= separated; // remove separated nodes from this partition
}

// check if two nodes in partition can separate
bool can_separate(GRAPH_T *graph, GRAPH_T a, GRAPH_T b, GRAPH_T *partitions,
                U_INT32 num_partitions)
{
#ifdef DEBUG
  printf("\n Check if unique");
#endif
  U_INT32 p=0;
  GRAPH_T mask;

  mask=ONE8<<(LEVEL-1); // position mask under first potential child

  //foreach child of a, b has a child in the same partition
  while(mask)
  {
    //if child in a exists, and it doesn't exists in b
    if((graph[a] & mask) && !(graph[b] & mask))
    { //child exists in a but not b
      for(p=0; p<num_partitions; ++p)
      {
        // if a's child is in partition[p]
        if(partitions[p] & mask)
        { // are any child of b in partition[p]?
          if(partitions[p] & graph[b])
            break;
          return true;
        }
      }
    }
    mask >>= 1;
  }
  return false;
}

// check for uniqueness by determining if this graph has duplicates
// through row swapping. Only count if this decoration
// has the highest value of all possible duplicates
bool is_unique(GRAPH_T *graph, U_INT64 orig)
{
#if LEVEL > 3
  GRAPH_T equiv_graph[LEVEL];
  GRAPH_T indexes[LEVEL]; // array holds current permute of middle row indexes
  U_INT32 i;              // interators
  U_INT64 perm;

  // create copy of graph to determine equivalent permutations
  memcpy(equiv_graph,graph,LEVEL*sizeof(GRAPH_T));

  // create initial permutation of {1,2,3,...,LEVEL-2}
  for (i=0; i<LEVEL; ++i)
    indexes[i]=i;

  while ((perm = next_permutation(equiv_graph, indexes)))
    if (perm > orig)
      return false;
#ifdef DEBUG
  if (perm == orig) printf("Hey! this equals that one! %llu\n",orig);
#endif
#endif
  return true;
}

// find next permutation of middle (LEVEL-2) rows from previous
// Use lexicographical ordering of an array containing the
// indexes for the main row (their original decoration)
// Every time we swap to create the next lexicographical permutation,
// we also swap the corresponding rows and columns in the graph
//
// Adapted from Kenneth Rosen's Discrete Mathematics and its Applications,
//  6th edition, p.384, Algorithm 1: Generating the Next Permutation
//  in Lexicographic Order
U_INT64 next_permutation(GRAPH_T *graph, GRAPH_T *indexes)
{
  GRAPH_T k,j,r,s;
  U_INT64 ans;

  k = LEVEL-3;
  while (indexes[k] > indexes[k+1]) k--;
  if (k == 0) return(0);
  else
  {
    j = LEVEL-2;
    while (indexes[k] > indexes[j]) j--;
    swap_rowcol(graph,indexes,j,k);
    r = LEVEL-2; s = k+1;
    while (r>s)
    {
      swap_rowcol(graph,indexes,r,s);
      r--; s++;
    }
  }
#ifdef DEBUG
  printf("\nNext permutation: ");
  for (j = 0; j<LEVEL; ++j)
    printf("%i ",indexes[j]);
  printf("\n");
#endif

  // compare values
  ans = graph2int(graph);
#ifdef DEBUG
  printf("Permutation is %llu\n",ans);
#endif
  return ans;
}

// swap a single row and column, and the indexes that refer to them
// for use with next_permutation()
void swap_rowcol(GRAPH_T *graph, GRAPH_T *indexes, GRAPH_T x, GRAPH_T y)
{
  GRAPH_T b1, b2, r;

  // swap indexes
  indexes[x] ^= indexes[y];
  indexes[y] ^= indexes[x];
  indexes[x] ^= indexes[y];

  // swap rows first
  graph[x] ^= graph[y];
  graph[y] ^= graph[x];
  graph[x] ^= graph[y];

  // swap columns x,y for i=0, i<LEVEL
  b1=ONE8<<(LEVEL-1-x);
  b2=ONE8<<(LEVEL-1-y);

  for (r=0; r<LEVEL; ++r)
  {
    // swap bits only if they are different
    // if they are different, the swap doesn't
    // need to know what they are: just that they
    // need to be flipped
    if(((graph[r] & b1)?1:0) != ((graph[r] & b2)?1:0))
      graph[r] ^= (b1 | b2);
  }
}

// create integer from graph
U_INT64 graph2int(GRAPH_T *graph)
{
  U_INT64 value=0;
  U_INT32 i;
  for (i=0; i<LEVEL; ++i)
  {
    value <<= LEVEL;
    value |= graph[i];
  }
  return value;
}
