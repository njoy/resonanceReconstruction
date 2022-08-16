#ifndef NJOY_R2_RMATRIX_LEGACY_RESOLVED_RESONANCE
#define NJOY_R2_RMATRIX_LEGACY_RESOLVED_RESONANCE

// system includes

// other includes
#include "resonanceReconstruction/Quantity.hpp"

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {
namespace legacy {
namespace resolved {

/**
 *  @class
 *  @brief Resolved SLBW or MLBW resonance parameters for a given energy
 *
 *  The resolved resonance class contains the width for elastic,
 *  capture, fission and total (given in eV) at a given energy.
 *
 *  The SLBW and MLBW ENDF resonance widths are defined as GT, GN, GG and GF.
 *  The energy dependent widths used in the resonance reconstruction formulas
 *  for these reactions are related to the stored widths as:
 *    gamma( E ) = P( E ) G( Er ) / P( Er )
 *  where E is the incidemt neutron energy and Er is the resonance energy.
 *
 *  Since the penetrability for fission and capture (and total which is only
 *  used to calculate the competitive width) is equal to 1, only the elastic
 *  width is energy dependent.
 */
class Resonance {

  /* fields */
  Energy energy_;

  Width elastic_;
  Width capture_;
  Width fission_;
  Width competition_;

  Width elastic_to_penetrability_;
  Width competition_to_penetrability_;
  double penetrability_;
  double shiftfactor_;

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/legacy/resolved/Resonance/src/ctor.hpp"

  /**
   *  @brief Return the resonance energy (in eV)
   */
  const Energy& energy() const { return this->energy_; }

  /**
   *  @brief Return the elastic width (in eV)
   */
  const Width& elastic() const { return this->elastic_; }

  /**
   *  @brief Return the capture width (in eV)
   */
  const Width& capture() const { return this->capture_; }

  /**
   *  @brief Return the fission width (in eV)
   */
  const Width& fission() const { return this->fission_; }

  /**
   *  @brief Return whether or not the resonance has fission or not
   */
  bool hasFission() const { return this->fission_.value != 0.0; }

  /**
   *  @brief Return the competitive width (in eV)
   */
  const Width& competition() const { return this->competition_; }

  /**
   *  @brief Return the total width (in eV)
   */
  Width total() const {

    return this->elastic() + this->capture()
           + this->fission() + this->competition();
  }

  /**
   *  @brief Return the neutron penetrability at the resonance energy
   */
  double penetrability() const { return this->penetrability_; }

  /**
   *  @brief Return the neutron shift factor at the resonance energy
   */
  double shiftfactor() const { return this->shiftfactor_; }

  /**
   *  @brief Return the elastic width (in eV)
   */
  Width elastic( double penetrability ) const {

    return penetrability * this->elastic_to_penetrability_;
  }

  /**
   *  @brief Return the competitive width (in eV)
   */
  Width competition( double penetrability ) const {

    return penetrability * this->competition_to_penetrability_;
  }

  /**
   *  @brief Return the total width (in eV)
   */
  Width total( double P, double Q ) const {

    return this->elastic( P ) + this->capture()
           + this->fission() + this->competition( Q );
  }

  /**
   *  @brief Return the primed resonance energy
   *
   *  The primed resonance energy is defined as:
   *    Eprime = Er + ( S(Er) - S(E) ) / ( 2 * P(Er) ) * GN
   */
  Energy energyPrime( double shiftfactor ) const {

    return this->energy()
           + ( this->shiftfactor_ - shiftfactor ) / 2.
             * this->elastic_to_penetrability_;
  }
};

} // resolved namespace
} // legacy namespace
} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
