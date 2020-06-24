/**
 *  @brief Return whether or not the energy is below the threshold for this
 *         channel
 *
 *  The incident energy is below the threshold energy for the channel if
 *      energy * ratio + q < 0.0
 *  where energy is the incident energy, ratio is the mass ratio M / ( m + M )
 *  for the incident particle pair and q is the Q value for this channel.
 *
 *  @param[in] energy   the energy to be tested
 */
bool belowThreshold( const Energy& energy ) {

  const Energy value = this->incidentParticlePair().massRatio() * energy;
  return ( value + this->Q() ) < 0.0 * electronVolt;
}
