/**
 *  @brief Evaluate the cross sections at the given energy
 *
 *  @param[in] energy       the incident energy
 *  @param[in,out] result   a map containing the accumulated cross sections
 */
void evaluate( const Energy& energy,
               std::map< ReactionID, CrossSection >& result ) {

  // data we need: k, P, phi, rho, g_J
  const auto channel = this->incidentChannel()
  const auto incident = channel.particlePair().pairID();
  const auto waveNumber = channel.waveNumber( energy );
  const auto penetrability = channel.penetrability( energy );
  const auto phaseShift = channel.phaseShift( energy );
  const auto ratio = waveNumber * channel.radii().penetrabilityRadius( energy );
  const auto spinFactor = channel.statisticalSpinFactor();
  const auto sinphi = std::sin( phaseShift );
  const auto sin2phi = sinphi * sinphi;

  // the 2 * pi2 / k2 factor
  const CrossSection factor =  2. * pi * pi / ( waveNumber * waveNumber );

  // interpolate on the resonance parameters at this energy and get the level
  // spacing and widths
  const auto parameters = this->resonanceTable()( energy );
  const auto spacing = parameters.levelSpacing();
  const Degrees degrees = this->resonanceTable().degreesOfFreedom();
  const Widths widths{ parameters.elastic() * sqrt( energy )
                         * degrees.elastic * penetrability / ratio ,
                       parameters.capture(),
                       parameters.fission(),
                       parameters.competition() };

  // calculate the fluctuation integrals
  const FluctuationIntegrals integrals =
    calculateFluctuationIntegrals( widths, degrees );

  // calculate the resulting cross sections
  result[ ReactionID( incident, incident ) ] +=
           factor * ( spinFactor / spacing *
             ( widths.elastic * widths.elastic * integrals.elastic
               - 2. * widths.elastic * sin2phi ) );
  result[ ReactionID( incident.symbol() + "->capture" ) ] +=
           factor * spinFactor / spacing *
             ( widths.elastic * widths.capture * integrals.capture );
  if ( widths.hasFission() ) {

    result[ ReactionID( incident.symbol() + "->fission" ) ] +=
             factor * spinFactor / spacing *
               ( widths.elastic * widths.fission * integrals.fission );
  }
}
