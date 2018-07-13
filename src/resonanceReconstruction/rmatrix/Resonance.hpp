/**
 *  @class
 *  @brief Resonance parameters in the Reich-Moore approximation with a single
 *         eliminated capture channel
 */
class Resonance {

  /* fields */
  Energy energy_;
  std::vector< ReducedWidth > widths_;
  ReducedWidth eliminated_;

  /* auxiliary functions */

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/Resonance/src/ctor.hpp"

  /**
   *  @brief Return the resonance energy (in eV)
   */
  auto energy() const {

    return this->energy_;
  }

  /**
   *  @brief Return the eliminated capture width (in sqrt(eV))
   */
  auto eliminatedWidth() const {

    return this->eliminated_;
  }

  /**
   *  @brief Return the reduced widths for this resonance (in sqrt(eV))
   */
  auto widths() const {

    return ranges::view::all( this->widths_ );
  }

  #include "resonanceReconstruction/rmatrix/Resonance/src/rmatrix.hpp"
};
