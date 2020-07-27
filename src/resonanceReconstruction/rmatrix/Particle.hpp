/**
 *  @class
 *  @brief Particle information
 *
 *  The Particle class contains specific information for a particle as used
 *  during resonance reconstruction. The Particle has an atomic mass, an 
 *  electrical charge, a spin and a parity (either + or -).
 *
 *  These variables are used to calculate quantities like the wave number k for
 *  an incident reaction channel.
 */
class Particle {

  /* fields */
  ParticleID id_;
  AtomicMass mass_;
  ElectricalCharge charge_;
  Spin spin_;
  Parity parity_;

  /* auxiliary functions */
  #include "resonanceReconstruction/rmatrix/Particle/src/verifyNotNegative.hpp"

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/Particle/src/ctor.hpp"

  /**
   *  @brief Return the particle identifier
   */
  const ParticleID& particleID() const { return this->id_; }

  /**
   *  @brief Return the atomic mass of the particle
   */
  const AtomicMass& mass() const { return this->mass_; }

  /**
   *  @brief Return the electrical charge of the particle
   */
  const ElectricalCharge& charge() const { return this->charge_; }

  /**
   *  @brief Return the spin of the particle
   */
  const Spin& spin() const { return this->spin_; }

  /**
   *  @brief Return the particle parity
   */
  const Parity& parity() const { return this->parity_; }
};
