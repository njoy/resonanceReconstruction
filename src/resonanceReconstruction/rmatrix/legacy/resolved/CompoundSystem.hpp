#ifndef NJOY_R2_RMATRIX_LEGACY_RESOLVED_COMPOUNDSYSTEM
#define NJOY_R2_RMATRIX_LEGACY_RESOLVED_COMPOUNDSYSTEM

// system includes

// other includes
#include "range/v3/range/conversion.hpp"
#include "range/v3/action/sort.hpp"
#include "range/v3/action/unique.hpp"
#include "resonanceReconstruction/rmatrix/identifiers.hpp"
#include "resonanceReconstruction/rmatrix/options.hpp"
#include "resonanceReconstruction/rmatrix/Map.hpp"
#include "resonanceReconstruction/rmatrix/legacy/CompoundSystemBase.hpp"
#include "resonanceReconstruction/rmatrix/legacy/resolved/SpinGroup.hpp"

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {
namespace legacy {
namespace resolved {

/**
 *  @class
 *  @brief The compound system for legacy resolved SLBW and MLBW
 *
 *  This class contains the resolved resonance parameters for legacy ENDF
 *  data using SLBW and MLBW. It can be used to generate cross section values
 *  at any energy but it also provides the minimal energy grid that should be
 *  used (this grid is generated on request by the user).
 */
template < typename Formalism >
class CompoundSystem : protected CompoundSystemBase< SpinGroup< Formalism > > {

public:

  /* constructor */
  using CompoundSystemBase< SpinGroup< Formalism > >::CompoundSystemBase;

  /* methods */
  using CompoundSystemBase< SpinGroup< Formalism > >::spinGroups;
  using CompoundSystemBase< SpinGroup< Formalism > >::evaluate;
  using CompoundSystemBase< SpinGroup< Formalism > >::reactionIDs;

  #include "resonanceReconstruction/rmatrix/legacy/resolved/CompoundSystem/src/grid.hpp"
};

} // resolved namespace
} // legacy namespace
} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
