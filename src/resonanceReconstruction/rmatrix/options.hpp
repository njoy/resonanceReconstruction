#ifndef NJOY_R2_RMATRIX_OPTIONS
#define NJOY_R2_RMATRIX_OPTIONS

// system includes

// other includes

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {

  // particle and channel types
  struct Neutron {};
  struct Photon {};
  struct ChargedParticle {};
  struct Fission {};

  // R-matrix boundary condition options
  struct ShiftFactor {}; // shift factor eliminates the boundary condition
  struct Constant {};    // the boundary condition is constannt

  // R-matrix formalism options
  struct SingleLevelBreitWigner {};
  struct MultiLevelBreitWigner {};
  struct ReichMoore {};
  struct GeneralRMatrix {};

} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
