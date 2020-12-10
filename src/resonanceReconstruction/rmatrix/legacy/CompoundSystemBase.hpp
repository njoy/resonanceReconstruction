#ifndef NJOY_R2_RMATRIX_LEGACY_COMPOUNDSYSTEMBASE
#define NJOY_R2_RMATRIX_LEGACY_COMPOUNDSYSTEMBASE

// system includes
#include <vector>

// other includes
#include "range/v3/algorithm/for_each.hpp"
#include "range/v3/algorithm/count_if.hpp"
#include "range/v3/view/all.hpp"
#include "range/v3/view/transform.hpp"
#include "resonanceReconstruction/Quantity.hpp"
#include "resonanceReconstruction/rmatrix/Map.hpp"
#include "resonanceReconstruction/rmatrix/ReactionID.hpp"
#include "resonanceReconstruction/rmatrix/calculatePhaseShift.hpp"

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {
namespace legacy {

/**
 *  @class
 *  @brief The base interface for a legacy compound system
 *
 *  This class contains the resolved or unresolved resonance parameters
 *  for legacy ENDF data.
 */
template < typename SpinGroupType >
class CompoundSystemBase {

  /* fields */
  std::vector< SpinGroupType > groups_;
  unsigned int lmax_;

  /* auxiliary functions */
  #include "resonanceReconstruction/rmatrix/legacy/CompoundSystemBase/src/getLMax.hpp"
  #include "resonanceReconstruction/rmatrix/legacy/CompoundSystemBase/src/verifySpinGroups.hpp"

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/legacy/CompoundSystemBase/src/ctor.hpp"

  /**
   *  @brief Return the l,J data
   */
  auto spinGroups() const { return ranges::view::all( this->groups_ ); }

  #include "resonanceReconstruction/rmatrix/legacy/CompoundSystemBase/src/evaluate.hpp"
};

} // legacy namespace
} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
