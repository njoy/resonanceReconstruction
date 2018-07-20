/**
 *  @brief Return the value of eta for the particle pair at a given energy
 *
 *  The parameter eta is an energy dependent quantity defined as follows:
 *     eta = za * zb * mu / hbar / k
 *  in which za and zb are the electrical charge of the particles in the 
 *  particle pair, mu is the reduced mass of the particle pair, hbar is the 
 *  Planck constant and k is the wave number.
 *
 *  @param energy   the energy at which the permeability needs to be evaluated
 */
EtaParameter etaParameter( const Energy& energy ) const {

  const auto za = this->particle().charge();
  const auto zb = this->residual().charge();
  const auto mu = this->reducedMass();
  const auto k = this->waveNumber( energy );
  return ( za * zb * mu ) / ( hbar * k ) ;
}

