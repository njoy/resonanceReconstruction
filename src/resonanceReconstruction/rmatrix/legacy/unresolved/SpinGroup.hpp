#ifndef NJOY_R2_RMATRIX_LEGACY_UNRESOLVED_SPINGROUP
#define NJOY_R2_RMATRIX_LEGACY_UNRESOLVED_SPINGROUP

// system includes

// other includes
#include "resonanceReconstruction/Quantity.hpp"
#include "resonanceReconstruction/rmatrix/identifiers.hpp"
#include "resonanceReconstruction/rmatrix/Map.hpp"
#include "resonanceReconstruction/rmatrix/legacy/SpinGroupBase.hpp"
#include "resonanceReconstruction/rmatrix/legacy/unresolved/calculateFluctuationIntegrals.hpp"
#include "resonanceReconstruction/rmatrix/legacy/unresolved/ResonanceTable.hpp"

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {
namespace legacy {
namespace unresolved {

/**
 *  @class
 *  @brief Unresolved resonance data for a specific l,J spin group
 *
 *  This class contains the unresolved resonance parameters and the associated
 *  incident channel for legacy ENDF data.
 */
class SpinGroup : protected SpinGroupBase< unresolved::ResonanceTable > {

public:

  /* constructor */

  using SpinGroupBase::SpinGroupBase;

  /* methods */

  using SpinGroupBase::incidentChannel;
  using SpinGroupBase::orbitalAngularMomentum;
  using SpinGroupBase::totalAngularMomentum;
  using SpinGroupBase::resonanceTable;
  using SpinGroupBase::reactionIDs;
  using SpinGroupBase::hasFission;

  #include "resonanceReconstruction/rmatrix/legacy/unresolved/SpinGroup/src/evaluate.hpp"
};

} // unresolved namespace
} // legacy namespace
} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
