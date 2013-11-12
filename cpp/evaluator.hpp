#ifndef EVALUATOR_HPP
#define EVALUATOR_HPP

/**
 * @author pkambadu
 */

#include <utility>
#include <algorithm>
#include <numeric>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>

#include "utilities.hpp"

struct prng_seed_params_t {
  friend struct prng_seed_generator_t;

public:
  prng_seed_params_t(const int rand_seed, const int epoch, const int mpi_rank):
    rand_seed(rand_seed), epoch(epoch), mpi_rank(mpi_rank) { }

private:
  const int rand_seed;
  const int epoch;
  const int mpi_rank;
};

/**
 * Generates a PRNG seed. The seed is derived from:
 * + rand_seed: the command line parameter
 * + epoch: the SGD epoch number
 * + mpi rank
 * + task id
 */
struct prng_seed_generator_t {
  typedef prng_seed_params_t PrngSeedParams;

public:
  prng_seed_generator_t(const PrngSeedParams &seed_params, const int task_id):
    seed_params(seed_params), task_id(task_id) { }

  uint32_t operator()() const {
    std::vector<uint32_t> words;
    words.push_back(seed_params.rand_seed);
    words.push_back(seed_params.epoch);
    words.push_back(seed_params.mpi_rank);
    words.push_back(task_id);
    return std::accumulate(words.begin(), words.end(), (uint32_t)0, mix);
  }

private:
  static uint32_t mix(uint32_t x, uint32_t y) {
    return boost::random::mt19937(x ^ y)();
  }

private:
  const PrngSeedParams &seed_params;
  const int task_id;
};

/**
 * Generates row numbers at random from a range according to a Solver-specific
 * distribution.
 */
template<typename GraphType, typename Solver>
struct row_sampling_distribution_t {
private:
  typedef typename GraphType::VertexIndexType Vtx;
  typedef
  row_sampling_distribution_t<GraphType, Solver>
  RowSamplingDistribution;

public:
  static RowSamplingDistribution* make(const GraphType &A,
                                       const Vtx v0,
                                       const Vtx v1) {
    int n = A.get_num_vertices();
    std::vector<double> p(n);
    Solver::probabilities(A, p);
    return new RowSamplingDistribution(p, v0, v1);
  }

private:
  row_sampling_distribution_t(const std::vector<double> &p,
                              const Vtx v0,
                              const Vtx v1):
    dist(p.begin() + v0, p.begin() + v1), v0(v0) { }

public:
  template<typename Prng>
  int operator()(Prng &prng) const {
    return v0 + dist(prng);
  }

private:
  const boost::random::discrete_distribution<> dist;
  const Vtx v0;
};

/**
 * Improves the approximate solution x to Ax=y by producing a sequence of random
 * rows and making one small improvement step sequentially for each row.
 *
 * Template parameters:
 * GraphType: type for A
 * RandomAccessContainer: vector type for x and y
 * Solver: takes a row index and makes the corresponding improvement step
 *
 * Produces the improved x as a side effect.
 */
template <typename GraphType,
          typename RandomAccessContainer,   
          typename Solver>
struct evaluator_t {
private:
  typedef prng_seed_generator_t PrngSeedGenerator;

public:
  typedef prng_seed_params_t PrngSeedParams;
  typedef row_sampling_distribution_t<GraphType, Solver> RowSamplingDistribution;
  
private:
  /* passed in parameters */
  const GraphType *A;
  const RandomAccessContainer *y;
  RandomAccessContainer *x;
  const PrngSeedParams *seed_params;
  const RowSamplingDistribution *dist;

public:
  /* Constructor */
  evaluator_t (const GraphType *A,
               const RandomAccessContainer *y,
               RandomAccessContainer *x,
               const PrngSeedParams *seed_params,
               const RowSamplingDistribution *dist) :
    A(A), y(y), x(x), seed_params(seed_params), dist(dist) { }

  /**
   * @param[in] space The list of points in space that need to be iterated.
   */
  void operator() (const pfunc::space_1D& space) {
    PrngSeedGenerator seed_gen((*seed_params), space.begin());
    boost::random::mt19937 prng(seed_gen());
    /* Loop over each element in the range that is given to us */
    for (size_t i=space.begin(); i<space.end(); ++i) {
      const int row = (*dist)(prng);

      /* Create the update for X[row] using the provided solver */
      Solver::solve ((*A), (*y), (*x), row);
    }
  }

  /**
   * When requested for a split, simply create a new selector and return.
   */ 
  evaluator_t split () const {
    return evaluator_t (A, y, x, seed_params, dist);
  }
      
  /** 
   * While joining, simply select one of the two results --- this will not 
   * work if the result is INVALID!
   */
  void join (const evaluator_t& other) { }
};

#endif // EVALUATOR_HPP
