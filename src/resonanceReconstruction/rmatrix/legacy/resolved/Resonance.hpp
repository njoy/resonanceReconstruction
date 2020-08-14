/**
 *  @class
 *  @brief Resolved SLBW or MLBW resonance parameters for a given energy
 *
 *  The resolved resonance class contains the reduced width for elastic,
 *  capture, fission and total (given in sqrt(eV)) at a given energy.
 *
 *  The SLBW and MLBW ENDF resonance widths are defined as GT, GN, GG and GF.
 *  The energy dependent widths used in the resonance reconstruction formulas
 *  for these reactions are related to the stored widths as:
 *    gamma( E ) = P( E ) G( Er ) / P( Er )
 *  where E is the incidemt neutron energy and Er is the resonance energy.
 *
 *  Since the penetrability for fission and capture (and total which is only
 *  used to calculte the competitive width) is equal to 1, only the elastic
 *  width is energy dependent.
 *
 *  The values stored here for a given Er are GT, GN and GF for total, capture
 *  and fission. For elastic, the value stored here should be GN / P( Er ).
 */
class Resonance {

  /* fields */
  Energy energy_;

  ReducedWidth total_;
  ReducedWidth elastic_;
  ReducedWidth capture_;
  ReducedWidth fission_;

  ReducedWidth delta_;

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/legacy/resolved/Resonance/src/ctor.hpp"

  /**
   *  @brief Return the resonance energy (in eV)
   */
  const Energy& energy() const { return this->energy_; }

  /**
   *  @brief Return the elastic reduced width (in sqrt(eV))
   */
  const ReducedWidth& elastic() const { return this->elastic_; }

  /**
   *  @brief Return the elastic reduced width (in sqrt(eV))
   */
  ReducedWidth elastic( double penetrability ) const {

    return this->elastic() * penetrability;
  }

  /**
   *  @brief Return the capture reduced width (in sqrt(eV))
   */
  const ReducedWidth& capture() const { return this->capture_; }

  /**
   *  @brief Return the fission reduced width (in sqrt(eV))
   */
  const ReducedWidth& fission() const { return this->fission_; }

  /**
   *  @brief Return the total reduced width (in sqrt(eV))
   */
  const ReducedWidth& total() const { return this->total_; }

  /**
   *  @brief Return the competitive width (in sqrt(eV))
   */
  ReducedWidth competition( double penetrability ) const {

    return this->delta_ - this->elastic( penetrability );
  }
};
