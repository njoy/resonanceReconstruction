/**
 *  @class
 *  @brief Resolved resonance parameters
 *
 *  The Resonance class contains the reduced widths (given in sqrt(eV)) for a
 *  number of channels. This class is both compatible with the general R matrix
 *  theory  and with the Reich-Moore approximation with a single eliminated
 *  capture channel. When using the general R matrix theory, the eliminated
 *  width is zero.
 */
class Resonance : protected BaseResonance {

  /* fields */
  ReducedWidth eliminated_;

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/Resonance/src/ctor.hpp"

  using BaseResonance::energy;
  using BaseResonance::widths;

  /**
   *  @brief Return the eliminated capture width (in sqrt(eV))
   */
  const ReducedWidth& eliminatedWidth() const { return this->eliminated_; }

  #include "resonanceReconstruction/rmatrix/Resonance/src/rmatrix.hpp"
};
