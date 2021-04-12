#ifndef NJOY_R2_RMATRIX_COMPOUNDSYSTEM
#define NJOY_R2_RMATRIX_COMPOUNDSYSTEM

// system includes
#include <algorithm>
#include <complex>
#include <vector>

// other includes
#include "Log.hpp"
#include "range/v3/action/sort.hpp"
#include "range/v3/action/unique.hpp"
#include "range/v3/algorithm/count_if.hpp"
#include "range/v3/algorithm/for_each.hpp"
#include "range/v3/view/all.hpp"
#include "range/v3/view/transform.hpp"
#include "resonanceReconstruction/Quantity.hpp"
#include "resonanceReconstruction/rmatrix/Map.hpp"
#include "resonanceReconstruction/rmatrix/ReactionID.hpp"
#include "resonanceReconstruction/rmatrix/ReactionChannelID.hpp"
#include "resonanceReconstruction/rmatrix/Formalism.hpp"
#include "resonanceReconstruction/rmatrix/BoundaryOption.hpp"
#include "resonanceReconstruction/rmatrix/SpinGroup.hpp"

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {

/**
 *  @class
 *  @brief The compound nucleus system
 */
template < typename Formalism, typename BoundaryOption >
class CompoundSystem {

  /* fields */
  std::vector< SpinGroup< Formalism, BoundaryOption > > groups_;
  std::vector< ReactionID > reactions_;

  /* auxiliary functions */
  #include "resonanceReconstruction/rmatrix/CompoundSystem/src/makeReactionIDs.hpp"
  #include "resonanceReconstruction/rmatrix/CompoundSystem/src/makeSpinGroups.hpp"
  #include "resonanceReconstruction/rmatrix/CompoundSystem/src/verifySpinGroups.hpp"

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/CompoundSystem/src/ctor.hpp"

  auto spinGroups() const {

    return ranges::cpp20::views::all( this->groups_ );
  }

  /**
   *  @brief Return the reactions defined in the compound system
   */
  auto reactionIDs() const {

    return ranges::cpp20::views::all( this->reactions_ );
  }

  //#include "resonanceReconstruction/rmatrix/CompoundSystem/src/switchIncidentPair.hpp"
  #include "resonanceReconstruction/rmatrix/CompoundSystem/src/evaluate.hpp"
  #include "resonanceReconstruction/rmatrix/CompoundSystem/src/evaluateTMatrix.hpp"
  #include "resonanceReconstruction/rmatrix/CompoundSystem/src/grid.hpp"
};

} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
