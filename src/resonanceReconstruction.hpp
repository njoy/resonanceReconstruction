#ifndef NJOY_RESONANCE_RECONSTRUCTION
#define NJOY_RESONANCE_RECONSTRUCTION

// system includes
#include <complex>

// physical quantities used in resonance reconstruction
#include "resonanceReconstruction/Quantity.hpp"

// ENDF components
#include "resonanceReconstruction/endf.hpp"

// the interpolation library
#include "interpolation.hpp"

// ranges
#include "range/v3/numeric/accumulate.hpp"
#include "range/v3/view/group_by.hpp"
#include "range/v3/view/transform.hpp"

// Eigen components
#include "Eigen/Dense"
#include "Eigen/LU"

namespace njoy {
namespace resonanceReconstruction {

using Matrix3x3 = Eigen::Matrix3cd;

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

#include "resonanceReconstruction/breitWigner.hpp"
#include "resonanceReconstruction/reichMoore.hpp"

}
}

#endif

// the rmatrix namespace and its components
#include "resonanceReconstruction/rmatrix.hpp"
