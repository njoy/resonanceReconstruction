/**
 *  @brief Return the wave number of the particle pair at a given energy
 *
 *  The wave number k is an energy dependent quantity defined as follows:
 *     hbar^2 k^2 = 2 * mu * ( energy * mu / ma + q )
 *  in which mu is the reduced mass of the particle pair, ma is the atomic
 *  mass of the first particle in the particle pair, q is the Q value
 *  associated to the particle pair and hbar is the Planck constant.
 *
 *  @param energy   the energy at which the wave number needs to be evaluated
 */
WaveNumber waveNumber( const Energy& energy ) const {

  const auto mu = this->reducedMass();
  const auto q = this->Q();
  return sqrt( 2. * mu * ( std::abs( energy * this->massratio_ + q ) ) ) / hbar;
}
