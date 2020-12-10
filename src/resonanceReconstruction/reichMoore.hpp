#ifndef NJOY_R2_REICHMOORE
#define NJOY_R2_REICHMOORE

// system includes

// other includes
#include "Eigen/Dense"
#include "Eigen/LU"
#include "resonanceReconstruction/breitWigner.hpp"

namespace njoy {
namespace resonanceReconstruction {

using Matrix3x3 = Eigen::Matrix3cd;

namespace reichMoore {

using breitWigner::CrossSection;

#include "resonanceReconstruction/reichMoore/Resonance.hpp"
#include "resonanceReconstruction/reichMoore/Lvalue.hpp"

struct Both{};
struct Neither{};
struct Channel{};
struct Scattering{};

#include "resonanceReconstruction/reichMoore/Inspector.hpp"
#include "resonanceReconstruction/reichMoore/Base.hpp"
#include "resonanceReconstruction/reichMoore/Type.hpp"
#include "resonanceReconstruction/reichMoore/Apply.hpp"

} // reichMoore namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
