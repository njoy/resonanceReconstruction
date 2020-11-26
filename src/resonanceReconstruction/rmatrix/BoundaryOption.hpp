#ifndef NJOY_R2_RMATRIX_BOUNDARYOPTION
#define NJOY_R2_RMATRIX_BOUNDARYOPTION

// system includes

// other includes

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {

  // R-matrix boundary condition options
  struct ShiftFactor {}; // shift factor eliminates the boundary condition
  struct Constant {};    // the boundary condition is constannt

} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
