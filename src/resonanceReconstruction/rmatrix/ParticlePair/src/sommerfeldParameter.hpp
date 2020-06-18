/**
 *  @brief Return the value of Sommerfeld parameter for the particle pair at a
 *         given energy
 *
 *  The Sommerfeld parameter eta is an energy dependent quantity defined as
 *  follows:
 *     eta = za * zb * mu / ( 4 * pi * epsilon0 * hbar^2 * k )
 *  in which za and zb are the electrical charge of the particles in the
 *  particle pair, mu is the reduced mass of the particle pair, hbar is the
 *  Planck constant, k is the wave number and epsilon0 is the vacuum
 *  permittivity.
 *
 *  It is a dimensionless parameter.
 *
 *  @param energy   the energy at which the Sommerfeld parameter needs to be evaluated
 */
double sommerfeldParameter( const Energy& energy ) const {

  using ElectronVoltMeter = decltype( electronVolt * meter );

  const auto za = this->particle().charge();
  const auto zb = this->residual().charge();
  const auto mu = this->reducedMass();
  const auto k = this->waveNumber( energy );
  return Quantity< ElectronVoltMeter >( za * zb / ( 4. * pi * epsilon0 ) )
           /  Quantity< ElectronVoltMeter >( ( hbar * hbar * k ) / mu );
}
