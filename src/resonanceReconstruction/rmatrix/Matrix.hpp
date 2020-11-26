#ifndef NJOY_R2_MATRIX
#define NJOY_R2_MATRIX

// system includes

// other includes
#include <Eigen/Core>
#include <Eigen/LU>

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {

  // matrix
  template < typename T > using Matrix = Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic >;
  template < typename T > using DiagonalMatrix = Eigen::DiagonalMatrix< T, Eigen::Dynamic >;

} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
