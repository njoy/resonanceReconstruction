/**
 *  @brief Return the wave number of the channel at a given energy
 *
 *  The wave number k is an energy dependent quantity defined as follows:
 *     hbar^2 k^2 = 2 * mu * ( energy * ratio + q )
 *  in which mu is the reduced mass of the channel's particle pair and ratio
 *  is the mass ratio M / ( m + M ) for the incident particle pair, q is the
 *  Q value associated to the transition of the incident particle pair to the
 *  channel's particle pair and hbar is the Planck constant.
 *
 *  @param energy   the energy at which the wave number needs to be evaluated
 */
WaveNumber waveNumber( const Energy& energy ) const {

  const auto mu = this->particlePair().reducedMass();
  const auto ratio = this->incidentParticlePair().massRatio();
  const auto q = this->Q();
  return sqrt( 2. * mu * ( std::abs( energy * ratio + q ) ) ) / hbar;
}
