// base classes
#include "resonanceReconstruction/rmatrix/legacy/SpinGroupBase.hpp"
#include "resonanceReconstruction/rmatrix/legacy/CompoundSystemBase.hpp"
#include "resonanceReconstruction/rmatrix/legacy/ResonanceTableBase.hpp"

// data collections for the legacy resolved and unresolved resonances
#include "resonanceReconstruction/rmatrix/legacy/Data.hpp"

#ifndef NJOY_R2_RMATRIX_LEGACY
#define NJOY_R2_RMATRIX_LEGACY

// system includes

// other includes

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {
namespace legacy {

  // legacy resolved and unresolved resonance reconstruction
  #include "resonanceReconstruction/rmatrix/legacy/resolved.hpp"
  #include "resonanceReconstruction/rmatrix/legacy/unresolved.hpp"

} // legacy namespace
} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
