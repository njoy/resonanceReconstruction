/**
 *  @class
 *  @brief Unresolved resonance parameters for a given energy
 *
 *  The unresolved resonance class contains the reduced width for elastic (given
 *  in sqrt(eV)) and the width for capture, fission and competition (given in
 *  eV) at a given energy. In addition to these widths, it also contains the
 *  average level spacing (given in eV).
 */
class Resonance {

  /* fields */
  Energy energy_;
  LevelSpacing spacing_;

  ReducedWidth elastic_;
  Width capture_;
  Width fission_;
  Width competition_;

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/legacy/unresolved/Resonance/src/ctor.hpp"

  /**
   *  @brief Return the resonance energy (in eV)
   */
  const Energy& energy() const { return this->energy_; }

  /**
   *  @brief Return the average level spacing (in eV)
   */
  const Energy& levelSpacing() const { return this->spacing_; }

  /**
   *  @brief Return the elastic reduced width (in sqrt(eV))
   */
  const ReducedWidth& elastic() const { return this->elastic_; }

  /**
   *  @brief Return the capture width (given in eV)
   */
  const Width& capture() const { return this->capture_; }

  /**
   *  @brief Return the fission width (given in eV)
   */
  const Width& fission() const { return this->fission_; }

  /**
   *  @brief Return the competitive width (given in eV)
   */
  const Width& competition() const { return this->competition_; }
};
