#ifndef NJOY_R2_RMATRIX_LEGACY_UNRESOLVED_COMPOUNDSYSTEM
#define NJOY_R2_RMATRIX_LEGACY_UNRESOLVED_COMPOUNDSYSTEM

// system includes

// other includes
#include "range/v3/action/sort.hpp"
#include "range/v3/action/unique.hpp"
#include "resonanceReconstruction/rmatrix/Map.hpp"
#include "resonanceReconstruction/rmatrix/Formalism.hpp"
#include "resonanceReconstruction/rmatrix/ReactionID.hpp"
#include "resonanceReconstruction/rmatrix/legacy/CompoundSystemBase.hpp"
#include "resonanceReconstruction/rmatrix/legacy/unresolved/SpinGroup.hpp"

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {
namespace legacy {
namespace unresolved {

/**
 *  @class
 *  @brief The legacy unresolved compound system
 *
 *  This class contains the unresolved resonance parameters for legacy ENDF
 *  data. It can be used to generate cross section values at any energy but it
 *  also provides the energy grid that should be used (this grid is generated
 *  on request by the user).
 */
class CompoundSystem : protected CompoundSystemBase< unresolved::SpinGroup > {

  int interpolation_;

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/legacy/unresolved/CompoundSystem/src/ctor.hpp"

  /* methods */
  using CompoundSystemBase::spinGroups;
  using CompoundSystemBase::evaluate;

  /**
   *  @brief Return the interpolation scheme to be applied
   */
  int interpolation() const { return this->interpolation_; }

  #include "resonanceReconstruction/rmatrix/legacy/unresolved/CompoundSystem/src/grid.hpp"
};

} // unresolved namespace
} // legacy namespace
} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
