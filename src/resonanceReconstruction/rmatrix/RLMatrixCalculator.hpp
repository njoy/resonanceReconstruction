#ifndef NJOY_R2_RMATRIX_RLMATRIXCALCULATOR
#define NJOY_R2_RMATRIX_RLMATRIXCALCULATOR

// system includes

// other includes
#include "range/v3/view/transform.hpp"
#include "resonanceReconstruction/Quantity.hpp"
#include "resonanceReconstruction/rmatrix/options.hpp"
#include "resonanceReconstruction/rmatrix/Matrix.hpp"
#include "resonanceReconstruction/rmatrix/Channel.hpp"
#include "resonanceReconstruction/rmatrix/ParticleChannel.hpp"
#include "resonanceReconstruction/rmatrix/ResonanceTable.hpp"
#include "resonanceReconstruction/rmatrix/LMatrixCalculator.hpp"

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {

template < typename Formalism, typename BoundaryOption >
class RLMatrixCalculator;

#include "resonanceReconstruction/rmatrix/RLMatrixCalculator/ReichMoore.hpp"
#include "resonanceReconstruction/rmatrix/RLMatrixCalculator/GeneralRMatrix.hpp"

} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
