#ifndef NJOY_R2_RMATRIX_LMATRIXCALCULATOR
#define NJOY_R2_RMATRIX_LMATRIXCALCULATOR

// system includes
#include <variant>
#include <complex>

// other includes
#include "range/v3/view/transform.hpp"
#include "range/v3/view/zip_with.hpp"
#include "resonanceReconstruction/quantities.hpp"
#include "resonanceReconstruction/Matrix.hpp"
#include "resonanceReconstruction/rmatrix/BoundaryOption.hpp"
#include "resonanceReconstruction/rmatrix/Channel.hpp"

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {

template < typename BoundaryOption > class LMatrixCalculator;

#include "resonanceReconstruction/rmatrix/LMatrixCalculator/Constant.hpp"
#include "resonanceReconstruction/rmatrix/LMatrixCalculator/ShiftFactor.hpp"

} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
