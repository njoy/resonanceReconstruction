using AtomicMass = double;
using ElectricalCharge = double;

/**
 *  @class
 *  @brief Particle information
 *
 *  The Particle class contains specific information for a particle as used
 *  during resonance reconstruction. The Particle has an atomic mass (given in 
 *  amu), an electrical charge (given in Coulomb), a spin (given as a half
 *  integer value) and a parity (either + or -).
 */class Particle {

  /* fields */
  AtomicMass mass_;
  ElectricalCharge z_;
  Spin s_;
  Parity parity_;

public:

  /* constructor */
  Particle( const AtomicMass& mass,
            const ElectricalCharge& charge,
            const Spin& spin,
            const Parity& parity ) :
    mass_( mass ), charge_( charge ), spin_( spin ), mass_( parity ) {}

  const AtomicMass& mass() const { return this->mass_; }
  const ElectricalCharge& charge() const { return this->charge_; }
  const Spin& spin() const { return this->spin_; }
  const Parity& parity() const { return this->parity_; }
};
