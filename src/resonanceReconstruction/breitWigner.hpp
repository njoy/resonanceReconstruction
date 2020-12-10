#ifndef NJOY_R2_BREITWIGNER
#define NJOY_R2_BREITWIGNER

// system includes
#include <complex>

// other includes
#include "range/v3/numeric/accumulate.hpp"
#include "range/v3/view/group_by.hpp"
#include "range/v3/view/transform.hpp"
#include "interpolation.hpp"
#include "resonanceReconstruction/Quantity.hpp"
#include "resonanceReconstruction/endf.hpp"

namespace njoy {
namespace resonanceReconstruction {

template< int i >
using Integer = std::integral_constant< int, i >;

#include "resonanceReconstruction/src/radius.hpp"
#include "resonanceReconstruction/src/channelRadius.hpp"
#include "resonanceReconstruction/src/root.hpp"
#include "resonanceReconstruction/src/neutronWaveNumber.hpp"
#include "resonanceReconstruction/src/penetrationShift.hpp"
#include "resonanceReconstruction/src/phaseShift.hpp"

#include "resonanceReconstruction/pack.hpp"
#include "resonanceReconstruction/EnergyRange.hpp"
#include "resonanceReconstruction/ZeroWidth.hpp"

namespace breitWigner {

#include "resonanceReconstruction/breitWigner/src/psiChi.hpp"
#include "resonanceReconstruction/breitWigner/CrossSection.hpp"

#include "resonanceReconstruction/breitWigner/resonance.hpp"
#include "resonanceReconstruction/breitWigner/lvalue.hpp"

#include "resonanceReconstruction/breitWigner/Type.hpp"

#include "resonanceReconstruction/breitWigner/singleLevel.hpp"
#include "resonanceReconstruction/breitWigner/multiLevel.hpp"

} // breitWigner namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
