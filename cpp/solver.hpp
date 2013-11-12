#ifndef SOLVER_HPP
#define SOLVER_HPP

/**
 * @author pkambadu
 */

#include "lapack_prototypes.hpp"
#include "utilities.hpp"

/**
 * Solves Ax=y for x using the randomized Jacobi method. The method works by
 * making a sequence of steps, where each step selects a random row i and
 * changes x_i so that element i of the residual y-Ax becomes 0.
 *
 * Convergence rate analysis is in:
 * Greedy and randomized versions of the multiplicative Schwarz method
 * Michael Griebel and Peter Oswald
 * LAA 437:1596--1610, http://dx.doi.org/10.1016/j.laa.2012.04.052
 */
template <typename GraphType,
          typename RandomAccessContainer>
struct jacobi_solver_t {
  public:
  /**
   * One step of randomized Jacobi. Changes x according to:
   *   x(row) <- 1/A(row, row) * (b(row) - A(row, :) * x).
   *
   * @param[in] A A sparse graph representation (M,N).
   * @param[in] y Vector of dimension (M)
   * @param[inout] x Vector of dimension (N)
   * @param[in] row Row of A to be used
   */
  static void solve (const GraphType& A,
                     const RandomAccessContainer& y,
                     RandomAccessContainer& x,
                     int row) {
    const int row_begin = A.begin (row);
    const int row_end = A.end (row);

    const double y_row = y[row];
    const double A_row_row = A.get_weight(A.get_target_index(row,row));

    double A_row_transpose_x = 0.0;
    for (int i=row_begin; i<row_end; ++i) {
      const int col = A.get_target (i);
      if (col != row) {
        const double val = A.get_weight (i);
        A_row_transpose_x += (val*x[col]);
      }
    }

    x[row] = (1.0/A_row_row) * (y_row - A_row_transpose_x);
  }

  /**
   * Produces a length-n vector with the sampling probabilities of the rows.
   * That vector is (1 / trace(A)) * diag(A).
   */
  static void probabilities(const GraphType &A, RandomAccessContainer &p) {
    A.diag(p);
    double trace = std::accumulate(p.begin(), p.end(), 0.0);
    std::transform(p.begin(), p.end(),
                   p.begin(),
                   std::bind2nd(std::divides<double>(), trace));
  }
};

#endif // SOLVER_HPP
