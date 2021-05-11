#ifndef NJOY_R2_RMATRIX_RESONANCE
#define NJOY_R2_RMATRIX_RESONANCE

// system includes
#include <complex>
#include <vector>

// other includes
#include "range/v3/view/cartesian_product.hpp"
#include "range/v3/view/chunk.hpp"
#include "range/v3/view/transform.hpp"
#include "resonanceReconstruction/Quantity.hpp"

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {

/**
 *  @class
 *  @brief Resolved resonance parameters for a single resonance
 *
 *  The Resonance class contains the reduced widths (given in sqrt(eV)) for a
 *  number of channels. This class is both compatible with the general R matrix
 *  theory  and with the Reich-Moore approximation with a single eliminated
 *  capture channel. When using the general R matrix theory, the eliminated
 *  width is zero.
 */
class Resonance {

  /* fields */
  Energy energy_;
  std::vector< ReducedWidth > widths_;
  ReducedWidth eliminated_;

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/Resonance/src/ctor.hpp"

  /**
   *  @brief Return the resonance energy (in eV)
   */
  const Energy& energy() const { return this->energy_; }

  /**
   *  @brief Return the eliminated capture width (in sqrt(eV))
   */
  const ReducedWidth& eliminatedWidth() const { return this->eliminated_; }

  /**
   *  @brief Return the reduced widths for this resonance (in sqrt(eV))
   */
  auto widths() const {

    return ranges::cpp20::views::all( this->widths_ );
  }

  #include "resonanceReconstruction/rmatrix/Resonance/src/rmatrix.hpp"
};

} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
