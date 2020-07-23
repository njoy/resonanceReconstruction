/**
 *  @class
 *  @brief Base class for resonance parameters
 */
class BaseResonance {

  /* fields */
  Energy energy_;
  std::vector< ReducedWidth > widths_;

protected:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/BaseResonance/src/ctor.hpp"

public:

  /**
   *  @brief Return the resonance energy (in eV)
   */
  const Energy& energy() const { return this->energy_; }

  /**
   *  @brief Return the reduced widths for this resonance (in sqrt(eV))
   */
  auto widths() const { return ranges::view::all( this->widths_ ); }
};
