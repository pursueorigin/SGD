#ifndef UTILITIES_HPP
#define UTILITIES_HPP

/**
 * @author pkambadu
 */

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <algorithm>
#include <iterator>
#include <ext/hash_map>
#include <ext/hash_set>

/* Definitions for all the parameters and constants that we will ever use */
static const int ROOT                    = 0;
static const int STOP_ON_OBJ_GAIN        = 1;

static const int INVALID_RANK            =  -1;
static const int INVALID_SIZE            =  -1;
static const int MAX_PATH_LENGTH         =1024;

static const int MAX_EPOCHS_INDEX        = 0;
static const int EPOCH_SIZE_INDEX        = 1;
static const int M_INDEX                 = 2;
static const int DEBUG_INDEX             = 3;
static const int NUM_THREADS_INDEX       = 4;
static const int RAND_SEED_INDEX         = 5;
static const int NUM_INT_PARAMETERS      = 6;
extern int int_params[NUM_INT_PARAMETERS];

static const int NUM_DBL_PARAMETERS      = 0;
extern double dbl_params[NUM_DBL_PARAMETERS];

static const int A_FILE_PATH_INDEX        = 0;
static const int Y_FILE_PATH_INDEX        = 1;
static const int NUM_CHR_PARAMETERS       = 2;
extern const char* chr_params[NUM_CHR_PARAMETERS];

/* These are defined for an MPI execution */
extern int mpi_rank;
extern int mpi_size;

namespace std {
  /* Specialize std::equal_to for strings; used for parameter passing */
  template <> 
  struct equal_to<const char*> {
    bool operator()(const char* one, const char* two) const {
      return (0==strcmp(one,two));
    }
  };

  /* Specialize std::equal_to for pairs; used for parameter passing */
  template <typename T1, typename T2> 
  struct equal_to<std::pair<T1,T2> >{
    bool operator()(const std::pair<T1,T2>& one, 
                    const std::pair<T1,T2>& two) const {
      return (one.first==two.first && one.second==two.second);
    }
  };
}

/* Print out a vector stored in a one-dimensional array */
template <typename RandomAccessIterator>
static inline void print_vector (RandomAccessIterator v,
                                 const int& M, 
                                 const char* prefix="Not named") {
  printf ("Vector %s=\n", prefix);
  for (int i=0; i<M; ++i) printf("%6.6lf ", v[i]); printf("\n");
}

/* Some constants that are used for LAPACK/BLAS calls */
#pragma GCC diagnostic ignored "-Wunused-variable"
static char TRANS = 't';
static char NO_TRANS = 'n';
static char UPPER = 'u';
static char LEFT = 'l';
static char NO_DIAG = 'n';

static int ONE_STEP = 1;

static double PLUS_ONE = 1.0;
static double MINUS_ONE = -1.0;
static double ZERO = 0.0;

#endif // UTILITIES_HPP
