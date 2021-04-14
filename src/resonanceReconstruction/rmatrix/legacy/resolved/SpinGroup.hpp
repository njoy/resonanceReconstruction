#ifndef NJOY_R2_RMATRIX_LEGACY_RESOLVED_SPINGROUP
#define NJOY_R2_RMATRIX_LEGACY_RESOLVED_SPINGROUP

// system includes
#include <vector>

// other includes
#include "range/v3/range/conversion.hpp"
#include "range/v3/action/sort.hpp"
#include "range/v3/action/unique.hpp"
#include "range/v3/view/all.hpp"
#include "range/v3/view/transform.hpp"
#include "range/v3/numeric/accumulate.hpp"
#include "range/v3/view/zip_with.hpp"
#include "resonanceReconstruction/Quantity.hpp"
#include "resonanceReconstruction/rmatrix/options.hpp"
#include "resonanceReconstruction/rmatrix/Map.hpp"
#include "resonanceReconstruction/rmatrix/ReactionID.hpp"
#include "resonanceReconstruction/rmatrix/legacy/Data.hpp"
#include "resonanceReconstruction/rmatrix/legacy/SpinGroupBase.hpp"
#include "resonanceReconstruction/rmatrix/legacy/resolved/ResonanceTable.hpp"

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {
namespace legacy {
namespace resolved {

template < typename Formalism > class SpinGroup;

#include "resonanceReconstruction/rmatrix/legacy/resolved/SpinGroup/SingleLevelBreitWigner.hpp"
#include "resonanceReconstruction/rmatrix/legacy/resolved/SpinGroup/MultiLevelBreitWigner.hpp"

} // resolved namespace
} // legacy namespace
} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
