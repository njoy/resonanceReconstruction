/**
 *  @class
 *  @brief Unresolved resonance parameters for a given energy
 *
 *  The unresolved resonance class contains the reduced widths (given in
 *  sqrt(eV)) for a number of channels at a given energy. In addition to these
 *  widths, it also contains the average level spacing (given in eV).
 */
class UnresolvedResonance : protected BaseResonance {

  /* fields */
  Energy spacing_;

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/UnresolvedResonance/src/ctor.hpp"

  using BaseResonance::energy;
  using BaseResonance::widths;

  /**
   *  @brief Return the eliminated capture width (in sqrt(eV))
   */
  const Energy& levelSpacing() const { return this->spacing_; }
};
