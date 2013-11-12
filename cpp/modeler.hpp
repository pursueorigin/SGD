#ifndef MODELER_HPP
#define MODELER_HPP

/**
 * @author pkambadu
 */

#include <vector>
#include <utility>

#include <boost/scoped_ptr.hpp>

#include <pfunc/pfunc.hpp>
#include <pfunc/utility.h>
#include <pfunc/space_1D.hpp>
#include <pfunc/parallel_reduce.hpp>

#include "utilities.hpp"
#include "partitioner.hpp"
#include "evaluator.hpp"
#include "omp_parallel_reduce.hpp"
#include "residual.hpp"
#include "output.hpp"

template <typename GraphType,
          typename RandomAccessContainer,
          typename Solver
#if USE_PFUNC
          ,typename PFuncGeneratorType
#endif
          >
struct modeler_t {
  typedef typename GraphType::VertexIndexType VertexIndexType;
  typedef typename GraphType::EdgeWeightType EdgeWeightType;
  typedef partitioner_t<VertexIndexType> partitioner_type;
  typedef typename partitioner_type::space_type space_type;

  /* Required for creating the evaluator */
  typedef evaluator_t<GraphType,
                      RandomAccessContainer,
                      Solver> Evaluator;
  typedef typename Evaluator::PrngSeedParams PrngSeedParams;
  typedef typename Evaluator::RowSamplingDistribution RowSamplingDistribution;
  typedef boost::scoped_ptr<RowSamplingDistribution> RowSamplingDistributionPtr;

#if USE_PFUNC
  /** Define the PFunc subtype */
  typedef typename PFuncGeneratorType::attribute attribute;
  typedef typename PFuncGeneratorType::task task;
  typedef typename PFuncGeneratorType::taskmgr taskmgr;
  typedef typename PFuncGeneratorType::functor functor;
  typedef pfunc::parallel_reduce <PFuncGeneratorType, 
                                  Evaluator, 
                                  space_type> ReduceType;
#endif


  private:
  /* passed as parameters */
  GraphType& A;
  const RandomAccessContainer& y;
  RandomAccessContainer& x;
  const int M;
  const int MAX_EPOCHS;
  const int EPOCH_SIZE;
  const int rand_seed;
  const int debug;

#if USE_PFUNC
  taskmgr* global_taskmgr;
#endif

  public:
  /* Constructor */
  modeler_t (GraphType& A,
             const RandomAccessContainer& y,
             RandomAccessContainer& x
#if USE_PFUNC
             ,taskmgr* global_taskmgr
#endif
             ) :
    A(A), y(y), x(x),
    M(int_params[M_INDEX]),
    MAX_EPOCHS(int_params[MAX_EPOCHS_INDEX]),
    EPOCH_SIZE(int_params[EPOCH_SIZE_INDEX]),
    rand_seed(int_params[RAND_SEED_INDEX]),
    debug(int_params[DEBUG_INDEX])
#if USE_PFUNC
    ,global_taskmgr(global_taskmgr)
#endif
  { }

  /**
   * This process basically synchronizes globally and makes 'x' consistent.
   * NOTE: We are assuming that EdgeWeightType is "double". This is a short
   * cut for now that needs to be removed. 
   * NOTE: We don't need to send the entire portion of 'x' to each and every
   * process. Just send those elements that are needed.
   */
  void global_sync (const space_type& my_row_space) {
    double* sendbuf = &(x[0]) + my_row_space.begin();
    int sendcount = my_row_space.end()-my_row_space.begin();
    double* recvbuf = &(x[0]);
    std::vector<int>recvcounts (mpi_size);
    std::vector<int>displs (1+mpi_size);

    /* Get the counts that each mpi_rank will send */
    partitioner_type::counts(0, M, mpi_size, recvcounts.begin());
    partitioner_type::intervals(0, M, mpi_size, displs.begin());

    MPI_Allgatherv (MPI_IN_PLACE, 0, MPI_DOUBLE,
                    recvbuf, &(recvcounts[0]), &(displs[0]), MPI_DOUBLE, 
                    MPI_COMM_WORLD);
  }

  /**
   * This is the main loop for the SGD iteration.
   *
   * Solves Ax=y for x by running an asynch solver on each thread on each MPI
   * rank. Runs a sequence of MAX_EPOCHS. Stops the computation after each
   * epoch, synchronizes x between all MPI ranks and then measures how good of a
   * solution it is. Each MPI rank makes a total of EPOCH_SIZE asynch steps in
   * each epoch.
   */
  void operator()() {
    /**
     * my_row_space: This is the range of random numbers to produce.
     * epoch_space: This is the number of rows to sample (shared by tasks).
     */
    space_type my_row_space 
      (partitioner_type::create (0, M, mpi_rank, mpi_size));
#if 0
    RowSamplingDistributionPtr pdist(RowSamplingDistribution::make
                                     (A,
                                      my_row_space.begin(),
                                      my_row_space.end()));
#else
    RowSamplingDistribution* pdist = 
      RowSamplingDistribution::make (A, my_row_space.begin(),
                                      my_row_space.end());
#endif

    output_formatter_t fmt;
    if (ROOT==mpi_rank && 2<debug) {
      fmt.print_hdr();
      fmt.print_progress(0, residual_norm(A, x, y));
    }

    /**
     * Iterate over the space [0..MAX_EPOCHS].
     */
    for (int epoch = 0; epoch<MAX_EPOCHS; ++epoch) {
      double epoch_time = micro_time();
      space_type epoch_space (partitioner_type::create (0, EPOCH_SIZE));
      const PrngSeedParams seed_params(rand_seed, epoch, mpi_rank);
      Evaluator evaluator (&A, &y, &x, &seed_params, pdist);
#if USE_PFUNC
      task root_task;
      attribute root_attribute (false /*nested*/, false /*grouped*/);
      ReduceType evaluate (epoch_space, evaluator, *global_taskmgr);
      pfunc::spawn (*global_taskmgr, root_task, root_attribute, evaluate);
      pfunc::wait (*global_taskmgr, root_task);
#else
      evaluator (epoch_space);
#endif

      /* Synchronize with global updates */
      global_sync (my_row_space);

      /* Print the results obtained in this epoch. */
      epoch_time = micro_time()-epoch_time;
      if (mpi_rank==ROOT && 2<debug) {
        printf ("Epoch %d: %lf (per sec)\n", 
                epoch,
                (mpi_size*(epoch_space.end()-epoch_space.begin()))/epoch_time);
      }

      if (ROOT==mpi_rank && 2<debug) 
        fmt.print_progress(epoch + 1, residual_norm(A, x, y));
    }
  }
};

#endif // MODELER_HPP
