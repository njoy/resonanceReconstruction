#ifndef NJOY_R2_MATRIX
#define NJOY_R2_MATRIX

// system includes

// other includes
#include <Eigen/Core>
#include <Eigen/LU>

namespace njoy {
namespace resonanceReconstruction {

  // matrix
  template < typename T > using Matrix = Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic >;
  template < typename T > using DiagonalMatrix = Eigen::DiagonalMatrix< T, Eigen::Dynamic >;

} // resonanceReconstruction namespace
} // njoy namespace

#endif
