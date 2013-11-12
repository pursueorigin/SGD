/**
 * @author pkambadu
 */

#include <iostream>
#include <functional>
#include <cstring>
#include <vector>
#include <list>
#include <utility>
#include <ext/hash_map>
#include <mpi.h>

#include <pfunc/pfunc.hpp>
#include <pfunc/utility.h>
#include <pfunc/space_1D.hpp>
#include <pfunc/parallel_reduce.hpp>

/****************************************/
/* Turn on/off PFunc --- ON by default */
#ifndef USE_PFUNC
# define USE_PFUNC 1
#endif
/****************************************/

#include "utilities.hpp"
#include "parser.hpp"
#include "matrix_market.hpp"
#include "reader.hpp"
#include "partitioner.hpp"
#include "weighted_csr_graph.hpp"
#include "modeler.hpp"
#include "solver.hpp"

/* These were declared as extern in utilities.hpp --- defining it here */
int int_params [NUM_INT_PARAMETERS];
double dbl_params [NUM_DBL_PARAMETERS];
const char* chr_params[NUM_CHR_PARAMETERS];
int mpi_rank = INVALID_RANK;
int mpi_size = INVALID_SIZE;

/* Type definitions required to get the weighted CSR graph going */
typedef int VertexIndexType;
typedef double EdgeWeightType;
typedef std::pair<VertexIndexType, EdgeWeightType> EdgeType;
typedef std::list<EdgeType> ForwardContainerType;
typedef std::vector<ForwardContainerType> AdjacencyListType;
typedef AdjacencyListType::iterator AdjacencyListTypeIterator;
typedef std::pair<VertexIndexType, VertexIndexType> VertexPairType;
typedef weighted_csr_graph<AdjacencyListType> GraphType;
typedef jacobi_solver_t<GraphType,std::vector<EdgeWeightType> > SolverType;

int main (int argc, char** argv) {
  /* Initialize MPI */
  MPI_Init (&argc, &argv);

  /* Figure out the rank and size */
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);

  /* MPI sends argc and argv everywhere --- parse everywhere */
  parse_parameters (argc,argv);

  /**
   * Now, we read 'A' --- for now, we are reading it on every single MPI 
   * place because we want to store the graph in a replicated manner. If 
   * you want to store the matrix in a partitioned manner, that can be 
   * done later by Alex.
   * TODO: Alex
   */
  AdjacencyListType adjacency_list;
                                                         
  std::pair<int,int> N_and_nnz = 
            read_matrix_market (chr_params[A_FILE_PATH_INDEX], 
                                true, /* has to be symmetric */
                                adjacency_list);
  size_t num_vertices = N_and_nnz.first;
  size_t num_edges = N_and_nnz.second;
  assert ((int)num_vertices == int_params[M_INDEX]);

  /* Store in Compressed Sparse Row format */
  GraphType A (adjacency_list, num_vertices, num_edges);
  A.sort ();

  /* de-allocate the adjacency_list */
  AdjacencyListTypeIterator iter = adjacency_list.begin ();
  AdjacencyListTypeIterator end = adjacency_list.end ();
  while (iter != end) (*iter++).clear();
  adjacency_list.clear ();

  /**
   * Now, we read in 'Y', which is assumed to be dense and stored in a plain
   * text file, one entry per line.
   * TODO: Alex
   */
  std::vector<EdgeWeightType> y(int_params[M_INDEX]);
  reader (chr_params[Y_FILE_PATH_INDEX], y.begin());

  /**
   * Create space for X.
   */
  std::vector<EdgeWeightType> x(int_params[M_INDEX],0.0);

  if (ROOT==mpi_rank && 2<int_params[DEBUG_INDEX]) {
    A.pretty_print();
    print_vector (y.begin(), int_params[M_INDEX], "Y");
  }

#if USE_PFUNC
  /**
   * Define the PFunc instance. Note that we HAVE TO USE PFUNC::USE_DEFAULT as
   * the type of the FUNCTOR so that we can use pfunc::parallel_reduce.
   */
  typedef
  pfunc::generator <pfunc::cilkS, /* Cilk-style scheduling */
                    pfunc::use_default, /* No task priorities needed */
                    pfunc::use_default /* any function type*/> generator_type;
  typedef generator_type::attribute attribute;
  typedef generator_type::task task;
  typedef generator_type::taskmgr taskmgr;

  /* Create an instance of PFunc if that is what is needed */
  taskmgr* global_taskmgr;
  const int n_queues = int_params [NUM_THREADS_INDEX];
  unsigned int* thds_per_q_arr = new unsigned int [n_queues];
  for (int i=0; i<n_queues; ++i) thds_per_q_arr [i] = ONE_STEP;
  global_taskmgr = new taskmgr (n_queues, thds_per_q_arr);
  delete [] thds_per_q_arr;
#endif

  /*************************************************************************/
  /*           Start the modeler and wait for good things                  */

  /* Concretize the modeler */
  typedef modeler_t<GraphType,  /* type for the sparse graph */
                    std::vector<EdgeWeightType>,  /* type for Y and X */
                    SolverType /* type for the solver */
#if USE_PFUNC
                    , generator_type  /* the generator type */
#endif
                    > my_modeler_t;
                      
  /* Create an instance of the modeler */
  my_modeler_t my_modeler (A,      /* data frame */
                           y,      /* regressor */
                           x       /* the output */
#if USE_PFUNC
                           ,global_taskmgr /* task manager for pfunc */
#endif
                           );

  /* Let the model compute */
  double time = micro_time ();
  my_modeler ();
  time = micro_time () - time;

  /* Print out any statistics that you want */
  if (ROOT==mpi_rank && 1<int_params[DEBUG_INDEX]) {
    printf ("Completed the task in %lf (secs)\n", time);
  }

#if USE_PFUNC
    delete global_taskmgr;
#endif

  /* Finalize MPI */
  MPI_Finalize ();

  return 0;
}
