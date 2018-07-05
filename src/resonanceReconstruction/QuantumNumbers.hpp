using OrbitalAngularMomentum = unsigned int;
using Spin = double;
using TotalAngularMomentum = double;
using Parity = short;

/**
 *  @class
 *  @brief The l,S,J,pi quantum numbers of a reaction channel
 *
 *  The QuantumNumbers class contains specific information for a particle as used
 *  during resonance reconstruction. The Particle has an atomic mass (given in 
 *  amu), an electrical charge (given in Coulomb), a spin (given as a half
 *  integer value) and a parity (either + or -).
 */
class QuantumNumbers {

  /* fields */
  OrbitalAngularMomentum l_;
  Spin s_;
  TotalAngularMomentum J_;
  Parity parity_;

public:

  /* constructor */
  QuantumNumbers( const OrbitalAngularMomentum& l,
                  const Spin& s,
                  const TotalAngularMomentum& J,
                  const Parity& parity ) :
    mass_( mass ), charge_( charge ), spin_( spin ), parity_( parity ) {}

  const OrbitalAngularMomentum& orbitalAngularMomentum() const { return this->l_; }
  const Spin& spin() const { return this->s_; }
  const TotalAngularMomentum& totalAngularMomentum() const { return this->J_; }
  const Parity& parity() const { return this->parity_; }
};
