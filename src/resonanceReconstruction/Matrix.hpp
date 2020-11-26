#ifndef NJOY_R2_MATRIX
#define NJOY_R2_MATRIX

// system includes

// other includes
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wextra"
#include <Eigen/Dense>
#pragma GCC diagnostic pop

namespace njoy {
namespace resonanceReconstruction {

  // matrix
  template < typename T > using Matrix = Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic >;
  template < typename T > using DiagonalMatrix = Eigen::DiagonalMatrix< T, Eigen::Dynamic >;

} // resonanceReconstruction namespace
} // njoy namespace

#endif
