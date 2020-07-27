/**
 *  @brief Return the value of Sommerfeld parameter for the channel at a
 *         given energy
 *
 *  The Sommerfeld parameter eta is an energy dependent quantity defined as
 *  follows:
 *     eta = z * Z * mu / ( 4 * pi * epsilon0 * hbar^2 * k )
 *  in which z and Z are the electrical charge of the particles in the
 *  particle pair, mu is the reduced mass of the particle pair, hbar is the
 *  Planck constant, k is the wave number and epsilon0 is the vacuum
 *  permittivity.
 *
 *  It is a dimensionless parameter.
 *
 *  @param energy   the energy at which the permeability needs to be evaluated
 */
double sommerfeldParameter( const Energy& energy ) const {

  using ElectronVoltMeter = decltype( electronVolt * meter );

  const auto zZ = this->particlePair().particle().charge() *
                  this->particlePair().residual().charge();
  if ( zZ.value == 0.0 ) {

    return 0.0;
  }
  else {

    const auto mu = this->particlePair().reducedMass();
    const auto k = this->waveNumber( energy );
    return Quantity< ElectronVoltMeter >( zZ / ( 4. * pi * epsilon0 ) )
           / Quantity< ElectronVoltMeter >( ( hbar * hbar * k ) / mu );
  }
}
