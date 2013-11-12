#ifndef RESIDUAL_HPP
#define RESIDUAL_HPP

#include <algorithm>
#include <numeric>

/**
 * Returns norm(y - Ax) / norm(y), where norm(.) is Euclidean norm.
 */
template<typename GraphType, typename RandomAccessContainer>
double residual_norm(const GraphType &A,
                     const RandomAccessContainer &x,
                     const RandomAccessContainer &y) {
  RandomAccessContainer r(A.get_num_vertices());
  std::transform(y.begin(), y.end(),
                 r.begin(),
                 std::bind1st(std::multiplies<double>(), -1.0));
  A.matvec(r, x);
  return std::sqrt(std::inner_product(r.begin(), r.end(),
                                      r.begin(),
                                      0.0)
                   /
                   std::inner_product(y.begin(), y.end(),
                                      y.begin(),
                                      0.0));
}

#endif // RESIDUAL_HPP
