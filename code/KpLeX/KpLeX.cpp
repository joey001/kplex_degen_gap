#define ALL
int KX;
#include <assert.h>
#include <bitset>
#include <iostream>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/times.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <algorithm>
using namespace std;

#define NONE -1
#define DELIMITER 0
#define PASSIVE 0
#define ACTIVE 1
#define MAX_NODE 80000000
#define max_expand_depth 100000
#define pop(stack) stack[--stack##_fill_pointer]
#define push(item, stack) stack[stack##_fill_pointer++] = item
#define ptr(stack) stack##_fill_pointer

#define CUR_KPX_SIZE KPX_Stack_fill_pointer
#define CURSOR Cursor_Stack[Cursor_Stack_fill_pointer - 1]

static int TIME_OUT, CUT_OFF = 0;
static double BEST_SOL_TIME;

static int NB_NODE, NB_NODE_O, NB_EDGE, NB_EDGE_O, MAX_KPX_SIZE, INIT_KPX_SIZE,
    HEUR_KPX_SIZE, MAX_VERTEX_NO;

#define CORE_NO Vertex_UB
static int Max_Degree = 0, UPPER_BOUND = 0;
static int Node_Degree[MAX_NODE];
static int *CNN;
static long long BRANCHING_COUNT = 0;

std::bitset<MAX_NODE> Node_State;
std::bitset<MAX_NODE> Node_State2;

static int **Node_Neibors;

static int Candidate_Stack_fill_pointer = 0;
static int Candidate_Stack[MAX_NODE * 2];
static int Vertex_UB[MAX_NODE * 2];
static int NBNN_Stack_fill_pointer = 0;
static char NBNN_Stack[MAX_NODE * 2];
static int  KPX_Stack_fill_pointer, Energy_Stack_fill_pointer;
static int *KPX_Stack, *MaxKPX_Stack, *Energy_Stack;
static int Cursor_Stack[max_expand_depth];
static int Cursor_Stack_fill_pointer = 0;

static int Rollback_Point;
static int Branching_Point;
static int *CORE_ORI;
static int NB_CANDIDATE = 0;

static int Extra_Node_Stack_fill_pointer = 0;
static int Extra_Node_Stack[1000000];

static int cut_ver = 0, total_cut_ver = 0;
static int cut_inc = 0, total_cut_inc = 0;
static int cut_iset = 0, total_cut_iset = 0;
static int cut_satz = 0, total_cut_satz = 0;

static int *Init_Adj_List;
static int BLOCK_COUNT = 0;
static int *BLOCK_LIST[50];
static double READ_TIME, INIT_TIME, SEARCH_TIME;
static double Dense0 = 0, Dense1 = 0, Dense2 = 0;

static double get_utime() {
  struct rusage utime;
  getrusage(RUSAGE_SELF, &utime);
  return (double)(utime.ru_utime.tv_sec + (double)utime.ru_utime.tv_usec / 1000000);
}

static int int_cmp_asc(const void *a, const void *b) {
  return *((int *)a) - *((int *)b);
}

static void printSol() {
  SEARCH_TIME = READ_TIME + INIT_TIME;
  printf("S ");
  sort(MaxKPX_Stack,MaxKPX_Stack+MAX_KPX_SIZE);
  for (int i = 0; i < MAX_KPX_SIZE; i++) {
    printf("%d ", MaxKPX_Stack[i]+1);
  }
  printf("\n");
}

static int is_adjacent(int node1, int node2) {
  int neibor, *neibors;
  if (node1 > node2) {
    neibors = Node_Neibors[node1];
  } else {
    neibors = Node_Neibors[node2];
    node2 = node1;
  }
  for (neibor = *neibors; neibor < node2 && neibor != NONE;
       neibor = *(++neibors))
    ;
  return neibor == node2;
}

static void allcoate_memory_for_adjacency_list(int nb_node, int nb_edge) {
  Init_Adj_List = (int *)malloc((2 * nb_edge + nb_node) * sizeof(int));
  BLOCK_COUNT = 1;
  BLOCK_LIST[BLOCK_COUNT - 1] = Init_Adj_List;
  Node_Neibors[1] = Init_Adj_List;
  for (int i = 2; i <= NB_NODE; i++) {
    Node_Neibors[i] = Node_Neibors[i - 1] + Node_Degree[i - 1] + 1;
  }
}

static int build_simple_graph_instance(char *input_file) {
  FILE *fp_in = fopen(input_file, "r");
  using ui = unsigned int;
  int m, n;
  ui t;
  fread(&t, sizeof(ui), 1, fp_in);
  fread(&n, sizeof(ui), 1, fp_in); // the number of vertices
  fread(&m, sizeof(ui), 1,
        fp_in); // the number of edges (twice the acutal number).
  fread(&Node_Degree[1], sizeof(ui), n, fp_in);
  NB_NODE = n;
  NB_EDGE = m;
  NB_NODE_O = NB_NODE;
  NB_EDGE_O = NB_EDGE;
  printf("R the graph size is %d\n", NB_NODE);
  Node_Neibors = (int **)malloc((NB_NODE + 1) * sizeof(int *));
  allcoate_memory_for_adjacency_list(NB_NODE, NB_EDGE);
  Max_Degree = 0;
  for (int i = 1; i <= NB_NODE; ++i) {
    fread(Node_Neibors[i], sizeof(ui), Node_Degree[i], fp_in);
    for (int j = 0; j < Node_Degree[i]; ++j)
      Node_Neibors[i][j]++;
    Node_Neibors[i][Node_Degree[i]] = NONE;
    if (Node_Degree[i] > Max_Degree)
      Max_Degree = Node_Degree[i];
  }
  UPPER_BOUND = Max_Degree + KX;
  Dense0 = ((float)NB_EDGE * 2 / NB_NODE / (NB_NODE - 1));
  printf("R Instance Information: #node=%d #edge=%d density=%9.8f\n", NB_NODE,
         NB_EDGE, Dense0);
  READ_TIME = get_utime();
  printf("R the reading time is %4.2lf \n", READ_TIME);
  return true;
}

static bool sort_by_degeneracy_ordering() {
  int *degree_counter, *where;
  int neibor, *neibors, pos, h, k;
  INIT_KPX_SIZE = 0;
  printf("I computing an initial k-plex...\n");

  where = Candidate_Stack + NB_NODE + 1;
  degree_counter = Vertex_UB + NB_NODE + 1;
  memset(degree_counter, 0, (Max_Degree + 1) * sizeof(int));

  for (int node = 1; node <= NB_NODE; node++) {
    Vertex_UB[node] = Node_Degree[node]; // core-number
    degree_counter[Node_Degree[node]]++;
  }

  int sum = 0; // prefix sum
  for (int i = 0; i <= Max_Degree; i++) {
    k = degree_counter[i];
    degree_counter[i] = sum;
    sum += k;
  }
  // order by degree
  for (int node = 1; node <= NB_NODE; node++) {
    Candidate_Stack[pos = degree_counter[Node_Degree[node]]++] = node;
    where[node] = pos;
  }

  for (int i = Max_Degree; i > 0; i--) {
    degree_counter[i] = degree_counter[i - 1];
  }
  degree_counter[0] = 0;

  Candidate_Stack[NB_NODE] = DELIMITER;
  ptr(Candidate_Stack) = NB_NODE + 1;

  int idx = 0;
  while (idx < NB_NODE) {
    int node = Candidate_Stack[idx];

    if (idx + 1 < NB_NODE &&
        Node_Degree[node] == Node_Degree[Candidate_Stack[idx + 1]]) {
      degree_counter[Node_Degree[node]] = idx + 1;
    }
    if (Node_Degree[node] > MAX_VERTEX_NO)
      MAX_VERTEX_NO = Node_Degree[node];
    if (Node_Degree[node] >=
        NB_NODE - idx - KX) { // when KX=1, the remaining is a clique.
      BEST_SOL_TIME = get_utime();
      MAX_KPX_SIZE = HEUR_KPX_SIZE = INIT_KPX_SIZE = NB_NODE - idx;
      printf("I the upper bound of k-plex %d ...\n"
             "I the initial %d-plex  %d...\n",
             UPPER_BOUND, KX, INIT_KPX_SIZE);
      MaxKPX_Stack = (int *)malloc((UPPER_BOUND + 1) *
                                   sizeof(int)); // UB is given during reading
      KPX_Stack = (int *)malloc((UPPER_BOUND + 1) * sizeof(int));
      Energy_Stack = (int *)malloc((UPPER_BOUND + 1) * sizeof(int));
      memcpy(MaxKPX_Stack, Candidate_Stack + idx, INIT_KPX_SIZE * sizeof(int));
      break;
    }
    neibors = Node_Neibors[node];
    for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
      if (where[neibor] > idx) {
        pos = where[neibor];
        h = degree_counter[Node_Degree[neibor]];

        k = Candidate_Stack[h];

        Candidate_Stack[h] = neibor;
        where[neibor] = h;

        Candidate_Stack[pos] = k;
        where[k] = pos;

        degree_counter[Node_Degree[neibor]]++;

        Node_Degree[neibor]--;
        if (Node_Degree[neibor] != Node_Degree[Candidate_Stack[h - 1]]) {
          degree_counter[Node_Degree[neibor]] = h;
        }
      }
    }
    idx++;
  }
  if (UPPER_BOUND == INIT_KPX_SIZE) {
    printf("I find the maximum %d-plex in initial phase!\n", KX);
    return true;
  }
  return false;
}

static void store_maximum_kpx(int node, bool silent) {
  push(node, KPX_Stack);
  MAX_KPX_SIZE++;
  BEST_SOL_TIME = get_utime();
  memcpy(MaxKPX_Stack, KPX_Stack, MAX_KPX_SIZE * sizeof(int));

  ptr(NBNN_Stack) = ptr(Candidate_Stack) = NB_NODE + 1;
  ptr(Cursor_Stack) = 1;
  ptr(KPX_Stack) = 0;
  ptr(Energy_Stack) = 0;
  Rollback_Point = 0;
  Vertex_UB[CURSOR] = MAX_KPX_SIZE;
  if (!silent) printf("C %4d |%7d |%14lld| %8.2lf\n", MAX_KPX_SIZE, CURSOR, BRANCHING_COUNT, BEST_SOL_TIME);
  total_cut_ver += cut_ver;
  cut_ver = 0;
  total_cut_inc += cut_inc;
  cut_inc = 0;
  total_cut_iset += cut_iset;
  cut_iset = 0;
  total_cut_satz += cut_satz;
  cut_satz = 0;
}

static inline void reset_state() {
  for (int i = CURSOR + 1, cn = Candidate_Stack[i]; cn != DELIMITER;
       cn = Candidate_Stack[++i]) {
    Node_State.reset(cn);
    Node_State2.reset(cn);
  }
}

static inline void reset_state1() {
  for (int i = 0; i < ptr(KPX_Stack); i++) {
    if (Energy_Stack[i] < 0) {
      Energy_Stack[i] = -Energy_Stack[i] - 1;
    }
  }
  for (int i = CURSOR + 1, cn = Candidate_Stack[i]; cn != DELIMITER;cn = Candidate_Stack[++i]) {
    Node_State.reset(cn);
    Node_State2.reset(cn);
  }
}

static inline void reset_state2() {
  for (int i = 0; i < ptr(KPX_Stack); i++) {
    if (Energy_Stack[i] < 0) {
      Energy_Stack[i] = -Energy_Stack[i] - 1;
    }
  }
}

static int produce_subgraph0() {
  int i = CURSOR, j = 0, neibor, max = 0, *neibors, *neibors2;

  int start = ptr(Candidate_Stack);

  int bnode = Candidate_Stack[CURSOR];

  assert(ptr(Candidate_Stack) == ptr(NBNN_Stack));

  for (int cn = Candidate_Stack[i = CURSOR + 1]; cn != DELIMITER; cn = Candidate_Stack[++i]) {
    Node_State.reset(cn);
    Node_State2.set(cn);
    // CNN[cn]=0;
  }

  // mark neighbors of current branching vertex,and set state
  neibors = Node_Neibors[bnode];
  for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
    if (Node_State2[neibor])
      Node_State.set(neibor);
  }

  NB_CANDIDATE = 0;
  for (int cn = Candidate_Stack[i = CURSOR + 1]; cn != DELIMITER; cn = Candidate_Stack[++i]) {
    assert(cn != bnode);
    if (Node_State[cn]) {
      push(cn, Candidate_Stack);
      push(NBNN_Stack[i], NBNN_Stack);
      NB_CANDIDATE++;
    } else if ((NBNN_Stack[CURSOR] < KX - 1) && (NBNN_Stack[i] < KX - 1)) {
      push(cn, Candidate_Stack);
      push(NBNN_Stack[i] + 1, NBNN_Stack);
      NB_CANDIDATE++;
    }
  }
  push(DELIMITER, Candidate_Stack);
  push(DELIMITER, NBNN_Stack);

  if (NB_CANDIDATE + CUR_KPX_SIZE <= MAX_KPX_SIZE) {
    reset_state();
    return false;
  }

  for (int cn = Candidate_Stack[i = CURSOR + 1]; cn != DELIMITER;cn = Candidate_Stack[++i]) {
    Node_State.reset(cn);
    Node_State2.reset(cn);
  }
  for (int cn = Candidate_Stack[j = start]; cn != DELIMITER;cn = Candidate_Stack[++j]) {
    Node_State2.set(cn);
  }

  int count = 0;
  for (i = 0; i < ptr(KPX_Stack) - 1; i++) {
    if(!is_adjacent(KPX_Stack[i], bnode)) {
      assert(Energy_Stack[i] < KX - 1);
      Energy_Stack[i] = -(Energy_Stack[i] + 1);
      count++;
    }
  }
  assert(count == NBNN_Stack[CURSOR]);
  assert(ptr(Candidate_Stack) == ptr(NBNN_Stack));

  for (i = 0; i < ptr(KPX_Stack) - 1; i++) {
    assert(Energy_Stack[i] < KX);
    if ((Energy_Stack[i] < 0 && (Energy_Stack[i] == 1 - KX))) {
      int nn = KPX_Stack[i];

      neibors = Node_Neibors[nn];
      for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
        if (Node_State2[neibor]) Node_State.set(neibor);
      }

      for (int cn = Candidate_Stack[j = start]; cn != DELIMITER; cn = Candidate_Stack[++j]) {
        assert(nn != cn || nn + cn != 0);
        if (cn > 0 && Node_State[cn] == 0) {
          Candidate_Stack[j] = -cn;
          NB_CANDIDATE--;
        }
      }

      neibors = Node_Neibors[nn];
      for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
        if (Node_State2[neibor]) Node_State.reset(neibor);
      }

      if (NB_CANDIDATE + CUR_KPX_SIZE <= MAX_KPX_SIZE) {
        for (int cn = Candidate_Stack[j = start]; cn != DELIMITER;cn = Candidate_Stack[++j]) {
          if (cn > 0) Node_State2.reset(cn);
          else Node_State2.reset(-cn);
        }
        for (int i = 0; i < ptr(KPX_Stack); i++) {
          if (Energy_Stack[i] < 0) {
            Energy_Stack[i] = -Energy_Stack[i] - 1;
          }
        }
        return false;
      }
    }
  }

  int _count = 0;
  i = j = start;
  for (int cn = Candidate_Stack[i]; cn != DELIMITER;cn = Candidate_Stack[++i]) {
    if (cn > 0) {
      Candidate_Stack[j] = cn;
      NBNN_Stack[j] = NBNN_Stack[i];
      _count++;
      j++;
    } else {
      Node_State2.reset(-cn);
    }
  }
  assert(_count == NB_CANDIDATE);
  Candidate_Stack[j] = DELIMITER;
  NBNN_Stack[j] = DELIMITER;
  ptr(Candidate_Stack) = j + 1;
  ptr(NBNN_Stack) = j + 1;

  if (NB_CANDIDATE + CUR_KPX_SIZE > MAX_KPX_SIZE) {
    return ptr(Candidate_Stack) - 1 - (MAX_KPX_SIZE - CUR_KPX_SIZE);
  } else {
    reset_state1();
    return false;
  }
}

static int cut_by_iteration_partition() {
  int ub = 0, k = ptr(Candidate_Stack) - 1;
  int lb = MAX_KPX_SIZE - CUR_KPX_SIZE;
  int min_k;
  if(KX==1) return k - (lb - ub);
  assert(Candidate_Stack[k] == DELIMITER);

#ifdef INVE
  for (int i = ptr(Clique_Stack) - 1; i >= 0; i--) {
#else
  for (int i = 0; i < ptr(KPX_Stack); i++) {
#endif
    int cnode = KPX_Stack[i];
    int ub1 = KX - 1 - Energy_Stack[i];
    assert(ub1 >= 0);
    if (ub1 == 0) continue;
    int j = k - 1, p = k, count = 0;
    for (int pnode = Candidate_Stack[j]; pnode != DELIMITER; pnode = Candidate_Stack[--j]) {
      if (is_adjacent(cnode, pnode) == false) {
        if (j < p - 1) {
          int temp = Candidate_Stack[p - 1];
          Candidate_Stack[p - 1] = pnode;
          Candidate_Stack[j] = temp;
          temp = NBNN_Stack[p - 1];
          NBNN_Stack[p - 1] = NBNN_Stack[j];
          NBNN_Stack[j] = temp;
        }
        p--;
        count++;
      }
    }
    assert(j < p);
    assert(Candidate_Stack[j] == DELIMITER);

    min_k = j;

    if (count < ub1)
      ub1 = count;

    if (ub + ub1 <= lb) {
      ub = ub + ub1;
      k = p;
      if (k == min_k + 1)
        break;
    } else {
      return k - (lb - ub) > min_k ? k - (lb - ub) : min_k + 1;
    }
  }
  if (k == min_k + 1) {
    return k;
  } else {
    assert(k > min_k + 1);
    int x = k - (lb - ub) > min_k ? k - (lb - ub) : min_k + 1;
    return x;
  }
}

static int cut_by_common_neibor(int start) {
  int i, j, neibor, max = 0, *neibors;
  int bnode = Candidate_Stack[CURSOR];
  if (CUR_KPX_SIZE >= 3) {
    for (int cn = Candidate_Stack[i = start]; cn != DELIMITER;cn = Candidate_Stack[++i]) {
      Node_State.reset(cn);
      Node_State2.reset(cn);
    }
    for (i = 0; i < ptr(KPX_Stack); i++) Energy_Stack[i]=abs(Energy_Stack[i]);
    return ptr(Candidate_Stack) - 1 - (MAX_KPX_SIZE - (CUR_KPX_SIZE));
  }

  neibors = Node_Neibors[bnode];
  for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
    if (Node_State2[neibor])
      Node_State.set(neibor);
  }

  for (int cn = Candidate_Stack[i = start]; cn != DELIMITER;cn = Candidate_Stack[++i]) {
    CNN[cn] = 0;
    neibors = Node_Neibors[cn];
    for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
      if (Node_State[neibor]) CNN[cn]++;
    }

    int flag = Node_State[cn] ? 0 : 1;
    int min_NN = MAX_KPX_SIZE - (CUR_KPX_SIZE + 1) + 1 - CNN[cn];

    if (NBNN_Stack[CURSOR] + NBNN_Stack[i] + min_NN + flag > 2 * KX - 2) {
      Node_State.reset(cn);
      Candidate_Stack[i] = -cn;
      NB_CANDIDATE--;
    }
  }

  for (int cn = Candidate_Stack[i = start]; cn != DELIMITER;cn = Candidate_Stack[++i]) {
      cn=abs(cn);
      Node_State.reset(cn);
      Node_State2.reset(cn);
  }

  if (NB_CANDIDATE + CUR_KPX_SIZE > MAX_KPX_SIZE) {
    int j = start, count = 0;
    float cutted = 0.0;
    for (int cn = Candidate_Stack[i = start]; cn != DELIMITER;
         cn = Candidate_Stack[++i]) {
      if (cn > 0) {
        Candidate_Stack[j] = cn;
        NBNN_Stack[j] = NBNN_Stack[i];
        count++;
        j++;
      } else {
        cutted++;
      }
    }
    Candidate_Stack[j] = DELIMITER;
    ptr(Candidate_Stack) = j + 1;
    NBNN_Stack[j] = DELIMITER;
    ptr(NBNN_Stack) = j + 1;

    assert(count == NB_CANDIDATE);

    for (i = 0; i < ptr(KPX_Stack); i++) {
      if (Energy_Stack[i] < 0) {
        Energy_Stack[i] = -Energy_Stack[i];
      }
    }
    return ptr(Candidate_Stack) - 1 - (MAX_KPX_SIZE - (CUR_KPX_SIZE));
  } else {
    for (int i = 0; i < ptr(KPX_Stack); i++) {
      if (Energy_Stack[i] < 0) {
        Energy_Stack[i] = -Energy_Stack[i] - 1;
      }
    }
    return false;
  }
}

static void init_for_search() {
  int i, node;
  int neibor, neibor2, *neibors, *neibors2;
  cut_ver = 0;
  cut_inc = 0;
  cut_iset = 0;
  cut_satz = 0;
  total_cut_ver = 0;
  total_cut_inc = 0;
  total_cut_iset = 0;
  total_cut_satz = 0;

  MAX_KPX_SIZE = INIT_KPX_SIZE;
  ptr(KPX_Stack) = 0;
  ptr(Energy_Stack) = 0;
  ptr(Cursor_Stack) = 0;
  push(NB_NODE - MAX_KPX_SIZE - 1, Cursor_Stack);
  Rollback_Point = 0;
  BRANCHING_COUNT = 0;

  assert(ptr(NBNN_Stack) == ptr(Candidate_Stack));

  for (i = 0; i < ptr(Candidate_Stack) - 1; i++) {
    node = Candidate_Stack[i];
    Vertex_UB[i] = Node_Degree[node] + KX;
    NBNN_Stack[i] = 0;
  }
  NBNN_Stack[NB_NODE] = DELIMITER;
  ptr(NBNN_Stack) = ptr(Candidate_Stack);
}

static void allocate_memory() {
  CNN = (int *)malloc((MAX_VERTEX_NO + 1) * sizeof(int));
}

static void search_maxKPX(int cutoff, bool silent) {
  init_for_search();
  MAX_KPX_SIZE=max(2*KX-2,MAX_KPX_SIZE);
  if (!silent) {
    printf("C  ---------------------------------------\n");
    printf("C  Size|   Index|   NB_Branches|   Time(s)\n");
  }
  int bnode;
  while (CURSOR > 0) {
    if (CUT_OFF > 0 && get_utime() > CUT_OFF) { TIME_OUT = true; break; }
    
    bnode = Candidate_Stack[--CURSOR];
    if (bnode == DELIMITER) {
      ptr(NBNN_Stack) = ptr(Candidate_Stack) = CURSOR + 1;
      ptr(Cursor_Stack)--;
      ptr(KPX_Stack)--;
      ptr(Energy_Stack)--;
      int pop_node = KPX_Stack[ptr(KPX_Stack)];
      assert(Candidate_Stack[CURSOR]==KPX_Stack[ptr(KPX_Stack)]);
      for (int i = 0; i < ptr(KPX_Stack); i++) {
        if (!is_adjacent(KPX_Stack[i], pop_node)) {
          Energy_Stack[i]--;
        }
      }
      // Vertex_UB[CURSOR]=MAX_KPX_SIZE-CUR_KPX_SIZE;
    } else {

      /* if(bnode<0) {
        bnode=-bnode;
        Candidate_Stack[CURSOR]=-Candidate_Stack[CURSOR];
      }
        */
      if (MAX_KPX_SIZE == CUR_KPX_SIZE) {
        store_maximum_kpx(bnode, silent);
      }
      // else if(0 && Vertex_UB[CURSOR] <= MAX_KPX_SIZE - CUR_KPX_SIZE) {
      //  cut_ver++;
      // }
      else {
        BRANCHING_COUNT++;
        Rollback_Point = ptr(Candidate_Stack);
        // push k-plex stack
        push(bnode, KPX_Stack);
        push(NBNN_Stack[CURSOR], Energy_Stack);

        if ((Branching_Point = produce_subgraph0()) == 0 ||
            (Branching_Point = cut_by_common_neibor(Rollback_Point)) == 0 ||
            (Branching_Point = cut_by_iteration_partition()) == 0) {
          ptr(NBNN_Stack) = ptr(Candidate_Stack) = Rollback_Point;
          ptr(KPX_Stack)--;
          ptr(Energy_Stack)--;
          continue;
        }
        push(Branching_Point, Cursor_Stack);
      }
    }
  }

  SEARCH_TIME = get_utime();
  if (!silent) {
    /*printf(
                          "C
       -----------------------------------------------------------------------------\n");
          printf("C %4d |%7d |%8d %10d %10d %10d|%14lld %8.2lf\n", MAX_KPX_SIZE,
       CURSOR,cut_ver,cut_inc, cut_iset, cut_satz,BRANCHING_COUNT,SEARCH_TIME);
          total_cut_ver += cut_ver;
          total_cut_inc += cut_inc;
          total_cut_iset += cut_iset;
          total_cut_satz += cut_satz;*/
    printf("C  ---------------------------------------\n");
    printf("C %4d |%7d |%14lld| %8.2lf\n", MAX_KPX_SIZE, CURSOR,BRANCHING_COUNT, SEARCH_TIME);
  }
}

static int *Adj_List;
#define ORI_CORE Node_Degree

static void free_block() {
  for (int i = 0; i < BLOCK_COUNT; i++)
    free(BLOCK_LIST[i]);
}

static void reduce_instance_with_unsupport_property() {
  int level = 0;
  int p1, p2, p3, i, node, *neibors1, *neibors2;

  for (int i = 1; i <= NB_NODE; i++) {
    assert(Candidate_Stack[NB_NODE - i] == i);
  }

  printf("R Reducing at level %3d", level);fflush(stdout);
  do {
    printf("\b\b\b%3d", ++level);fflush(stdout);
    for (p1 = Candidate_Stack[i = 0]; p1 != DELIMITER; p1 = Candidate_Stack[++i]) {
      if (p1 <= 0) continue;
      neibors1 = Node_Neibors[p1];
      for (p2 = *neibors1; p2 != p1; p2 = *(neibors1 += 2)) {
        *(neibors1 + 1) = 0;
      }
    }
    // list all triangles <p1,p2,p3>,computer common neibors for an edge <u,v>
    for (p1 = Candidate_Stack[i = 0]; p1 != DELIMITER; p1 = Candidate_Stack[++i]) {
      if (p1 <= 0) continue;
      assert(p1 == Candidate_Stack[NB_NODE - p1]);
      neibors1 = Node_Neibors[p1];
      for (p2 = *neibors1; p2 < p1; p2 = *(neibors1 += 2)) {
        if (p2 > 0) {
          if (Candidate_Stack[NB_NODE - p2] < 0)
            *(neibors1) = -p2;
          else
            Node_State.set(p2);
        } else if ((-p2) > p1) {
          break;
        }
      }
      neibors1 = Node_Neibors[p1];
      for (p2 = *neibors1; (p2 < p1) && (-p2 < p1); p2 = *(neibors1 += 2)) {
        if (p2 < 0) continue;
        assert(Node_State[p2]),assert(Node_Degree[p2] + KX > INIT_KPX_SIZE);
        neibors2 = Node_Neibors[p2];
        for (p3 = *neibors2; p3 < p2; p3 = *(neibors2 += 2)) {
          if (p3 > 0 && Node_State[p3]) {
            Node_State2.set(p3);
            *(neibors1 + 1) = *(neibors1 + 1) + 1;
            *(neibors2 + 1) = *(neibors2 + 1) + 1;
          } else if ((-p3) > p2) {
            break;
          }
        }
        neibors2 = Node_Neibors[p1];
        for (p3 = *neibors2; p3 < p2; p3 = *(neibors2 += 2)) {
          if (p3 > 0 && Node_State2[p3]) {
            *(neibors2 + 1) = *(neibors2 + 1) + 1;
            Node_State2.reset(p3);
          } else if ((-p3) > p2) {
            break;
          }
        }
      }

      neibors1 = Node_Neibors[p1];
      for (p2 = *neibors1; p2 < p1; p2 = *(neibors1 += 2)) {
        if (p2 > 0) {
          Node_State.reset(p2);
        } else if ((-p2) > p1) {
          break;
        }
      }
    }

    ptr(Extra_Node_Stack) = 0;
    // check all edges <p1,p2>, break <p1,p2> if comm(p1,p2)+2k<=LB
    // push vertices that can be removed
    for (p1 = Candidate_Stack[i = 0]; p1 != DELIMITER; p1 = Candidate_Stack[++i]) {
      if (p1 <= 0) continue;
      neibors1 = Node_Neibors[p1];
      for (p2 = *neibors1; (p2 < p1) && (-p2 < p1); p2 = *(neibors1 += 2)) {
        if (p2 > 0 && *(neibors1 + 1) + 2 * KX <= INIT_KPX_SIZE) {
          Node_Degree[p1]--;
          if (Node_Degree[p1] + KX == INIT_KPX_SIZE) push(p1, Extra_Node_Stack);
          *neibors1 = -p2;
          neibors2 = Node_Neibors[p2];
          for (p3 = *neibors2; p3 != p2; p3 = *(neibors2 += 2)) {
            if (p3 == p1) {
              *neibors2 = -p1;
              Node_Degree[p2]--;
              if (Node_Degree[p2] + KX == INIT_KPX_SIZE) push(p2, Extra_Node_Stack);
              break;
            }
          }
          assert(p3 == p1);
        }
      }
    }

    // removed vertices recursively
    push(DELIMITER, Extra_Node_Stack);
    for (int p1 = Extra_Node_Stack[i = 0]; p1 != DELIMITER; p1 = Extra_Node_Stack[++i]) {
      assert(p1 > 0);
      assert(Node_Degree[p1] <= INIT_KPX_SIZE - KX);
      assert(Candidate_Stack[NB_NODE - p1] == p1);

      Candidate_Stack[NB_NODE - p1] = -p1;
      neibors1 = Node_Neibors[p1];
      for (int p2 = *neibors1; p2 != p1; p2 = *(neibors1 += 2)) {
        if (p2 > 0 && Node_Degree[p2] + KX > INIT_KPX_SIZE) {
          Node_Degree[p2]--;
          if (Node_Degree[p2] + KX == INIT_KPX_SIZE) {
            Extra_Node_Stack[ptr(Extra_Node_Stack) - 1] = p2;
            push(DELIMITER, Extra_Node_Stack);
          }
        }
      }
    }
  }while(ptr(Extra_Node_Stack) > 1);
  printf("\n");
}

static int reduce_instance2() {
  int i = 0, j = 0, nb_edge = 0;
  int node, *neibors, *neibors2, *addr;
  MAX_VERTEX_NO = 0;
  for (int p1 = Candidate_Stack[i = 0]; p1 != DELIMITER; p1 = Candidate_Stack[++i]) {
    if (p1 > 0) {
      Candidate_Stack[j++] = p1;
      Node_State.set(p1);
      Node_Degree[p1] = 0;
    }
  }
  NB_NODE = j;

  if (NB_NODE <= INIT_KPX_SIZE) {
    printf("I prove the optimal in reduction phase!\n");
    return true;
  }

  ptr(NBNN_Stack) = ptr(Candidate_Stack) = j + 1;
  NBNN_Stack[j] = Candidate_Stack[j] = DELIMITER;

  NB_EDGE = 0;
  neibors2 = Adj_List;
  for (int p1 = Candidate_Stack[i = 0]; p1 != DELIMITER; p1 = Candidate_Stack[++i]) {
    neibors = Node_Neibors[p1];
    // neibors2 = neibors;
    for (node = *neibors; node != p1; node = *(neibors += 2)) {
      // if(node<0)node=-node;
      if (node > 0 && Node_State[node]) {
        Node_Degree[p1]++;
        *neibors2 = node;
        neibors2++;
        NB_EDGE++;
      }
    }
    (*neibors2) = NONE;
    neibors2++;
  }
  Node_Neibors[Candidate_Stack[0]] = Adj_List;
  for (int p0 = Candidate_Stack[0], p1 = Candidate_Stack[i = 1];p1 != DELIMITER; p0 = p1, p1 = Candidate_Stack[++i]) {
    Node_Neibors[p1] = Node_Neibors[p0] + Node_Degree[p0] + 1;
    Node_State.reset(p1);
  }
  NB_EDGE = NB_EDGE >> 1;
  MAX_VERTEX_NO = Candidate_Stack[0];
  Dense2=((float)2 * NB_EDGE / NB_NODE / (NB_NODE - 1));
  printf("I the L2-Reduced graph #node %d #edge %d #density %9.8f\n", j, NB_EDGE, Dense2);
  return false;
}

static int reduce_instance_with_degree() {
  int i = 0, NEW_NB_NODE = 0, node, *neibors, *new_neibors, *addr;
  // Node_Degree is core number now
  while (Node_Degree[Candidate_Stack[i]] + KX <= INIT_KPX_SIZE) {
    Node_State[Candidate_Stack[i]] = PASSIVE;
    i++;
  }
  while (i < NB_NODE) {
    node = Candidate_Stack[i];
    Candidate_Stack[NEW_NB_NODE++] = node;
    Node_State[node] = ACTIVE;
    i++;
  }
  NB_NODE = NEW_NB_NODE;
  if (NB_NODE <= INIT_KPX_SIZE) {
    printf("I find the optimal solution in level-1 reduction.\n");
    return true;
  }
  NBNN_Stack[NEW_NB_NODE] = Candidate_Stack[NEW_NB_NODE] = DELIMITER;
  ptr(NBNN_Stack) = ptr(Candidate_Stack) = NEW_NB_NODE + 1;

  CORE_ORI = (int *)malloc((NB_NODE + 1) * sizeof(int));
  // renumber (reverse for heuristic)
  for (i = 0; i < NB_NODE; i++) {
    CORE_ORI[NB_NODE - i] = Candidate_Stack[i];
    ORI_CORE[Candidate_Stack[i]] = NB_NODE - i;
    Candidate_Stack[i] = NB_NODE - i;
  }
  NB_EDGE = 0;
  for (i = NB_NODE; i > 0; i--) {
    neibors = Node_Neibors[CORE_ORI[i]];
    new_neibors = neibors;
    int new_neibor = 0;
    for (node = *neibors; node != NONE; node = *(++neibors)) {
      if (Node_State[node] == ACTIVE) {
        (*new_neibors) = ORI_CORE[node];
        new_neibors++;
        new_neibor++;
      }
    }
    (*new_neibors) = NONE;
    NB_EDGE += new_neibor;
    qsort(Node_Neibors[CORE_ORI[i]], new_neibor, sizeof(int), int_cmp_asc);
    assert(new_neibor + KX > INIT_KPX_SIZE);
  }
  NB_EDGE = NB_EDGE >> 1;
  Adj_List = (int *)malloc((4 * NB_EDGE + NB_NODE) * sizeof(int));
  addr = Adj_List;
  MAX_VERTEX_NO = 0;
  for (i = NB_NODE; i > 0; i--) {
    Node_Degree[i] = 0;
    Node_State.reset(i);
    Node_State2.reset(i);
    neibors = Node_Neibors[CORE_ORI[i]];
    for (node = *neibors; node != NONE; node = *(++neibors)) {
      *(addr++) = node;
      *(addr++) = 0;
      Node_Degree[i]++;
    }
    *(addr++) = i;
    if (Node_Degree[i] > MAX_VERTEX_NO)
      MAX_VERTEX_NO = Node_Degree[i];
  }
  free_block();
  Node_Neibors[NB_NODE] = Adj_List;
  for (i = NB_NODE - 1; i > 0; i--) {
    Node_Neibors[i] = Node_Neibors[i + 1] + 2 * Node_Degree[i + 1] + 1;
  }
  Dense1 = ((float)NB_EDGE * 2 / NB_NODE / (NB_NODE - 1));
  printf("I the L1-Reduced graph #node %d #edge %d #density %9.8f\n", NB_NODE,
         NB_EDGE, Dense1);
  return false;
}

static bool preprocess() {
  bool optimal = sort_by_degeneracy_ordering();
  if (!optimal)
    optimal = reduce_instance_with_degree();
  double begin = get_utime();
  if (!optimal)
    reduce_instance_with_unsupport_property(), optimal = reduce_instance2();
  double end = get_utime();
  printf("I the strong time is %4.2lf \n", end - begin);
  INIT_TIME = get_utime() - READ_TIME;
  printf("I the initial time is %4.2lf \n", INIT_TIME);
  return optimal;
}

static void check_and_printSol() {
  int count = 0, count1 = 0, node1;
  if (MAX_KPX_SIZE > INIT_KPX_SIZE) {
    for (int i = 0; i < MAX_KPX_SIZE; i++) {
      count = 0;
      count1 = 0;
      node1 = MaxKPX_Stack[i];
      for (int j = 0; j < MAX_KPX_SIZE; j++) {
        if (MaxKPX_Stack[j] == node1)
          count1++;
        else if (is_adjacent(node1, MaxKPX_Stack[j]) == false)
          count++;
      }
      if (count > KX - 1)
        std::cout << "count " << count << " " << KX - 1 << endl;
      assert(count1 == 1);
      assert(count <= KX - 1);
      // assert(count==Energy_Stack[i]);
    }
  }
  printf("S ");
  for (int i = 0; i < MAX_KPX_SIZE; i++) {
    if (MAX_KPX_SIZE > INIT_KPX_SIZE) {
      MaxKPX_Stack[i]=CORE_ORI[MaxKPX_Stack[i]]-1;
    }
  }
  sort(MaxKPX_Stack,MaxKPX_Stack+MAX_KPX_SIZE);
  for (int i = 0; i < MAX_KPX_SIZE; i++) {
      printf("%d ", MaxKPX_Stack[i]);
  }
  printf("\n");
}
int main(int argc, char *argv[]) {
  printf("\n-----------------------------------------------------------------------------------------\n");
  KX = atoi(argv[2]);
  CUT_OFF = 9999999;
  printf("# Solving %d-plex in %s\n", KX, argv[1]);
  build_simple_graph_instance(argv[1]);
  if (!preprocess()) {
    allocate_memory();
    search_maxKPX(0, 0);
    // check_and_printSol();
  } else
    printSol();
  printf(">>%s |V| %d |E| %d MaxKPX %d HEURKPX %d Tree %lld Read_Time %4.2lf "
         "Init_Time %4.2lf Search_Time %4.2lf Total %4.2lf \\\\\n",
         argv[1], NB_NODE_O, NB_EDGE_O, MAX_KPX_SIZE, HEUR_KPX_SIZE,
         BRANCHING_COUNT, READ_TIME, INIT_TIME,
         SEARCH_TIME - READ_TIME - INIT_TIME, SEARCH_TIME);
  printf("#File=%s\n#K=%d\n", argv[1], KX);
  printf("#nodes=%d\n#nedges=%d\n\n", NB_NODE_O, NB_EDGE_O);
  printf("#MaxKPX=%d\n#HEURKPX=%d\n#Tree=%lld\n\n", MAX_KPX_SIZE, HEUR_KPX_SIZE, BRANCHING_COUNT);
  printf("#timeout=%d\n#optimal=%d\n\n#Read_Time=%4.2lf\n#Init_Time=%4.2lf\n#Search_Time=%4.2lf\n#Total=%4.2lf\n\n", 1-TIME_OUT, CUT_OFF, READ_TIME, INIT_TIME,
         SEARCH_TIME - READ_TIME - INIT_TIME, SEARCH_TIME);
  printf("-----------------------------------------------------------------------------------------\n\n\n\n\n");
  return 0;
}
