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

  // R-matrix boundary condition option types
  struct ShiftFactor {}; // shift factor eliminates the boundary condition
  struct Constant {};    // the boundary condition is constannt

  // R-matrix boundary condition option enumerator

  enum class Boundary : short {

    ShiftFactor = 0, // shift factor eliminates the boundary condition
    Constant         // the boundary condition is constannt
  };

  // R-matrix formalism option types
  struct SingleLevelBreitWigner {};
  struct MultiLevelBreitWigner {};
  struct ReichMoore {};
  struct GeneralRMatrix {};

  // R-matrix formalism option enumerator

  enum class Formalism : short {

    SingleLevelBreitWigner = 0,
    MultiLevelBreitWigner,
    ReichMoore,
    GeneralRMatrix
  };

} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
