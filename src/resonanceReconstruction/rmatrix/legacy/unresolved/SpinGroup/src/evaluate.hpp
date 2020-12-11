/**
 *  @brief Evaluate the cross sections at the given energy
 *
 *  @param[in] energy       the incident energy
 *  @param[in,out] result   a map containing the accumulated cross sections
 */
void evaluate( const Energy& energy,
               Map< ReactionID, CrossSection >& result ) {

  // data we need: k, P, phi, rho, g_J
  decltype(auto) channel = this->incidentChannel();
  const auto waveNumber = channel.waveNumber( energy );
  const auto penetrability = channel.penetrability( energy );
  const auto phaseShift = channel.phaseShift( energy );
  const auto radius = channel.radii().penetrabilityRadius( energy );
  const auto spinFactor = channel.statisticalSpinFactor();
  const auto ratio = waveNumber * radius;
  const auto sinphi = std::sin( phaseShift );
  const auto sin2phi = sinphi * sinphi;

  // the 2 * pi2 / k2 factor
  const CrossSection factor = 2. * pi * pi / ( waveNumber * waveNumber );

  // interpolate on the resonance parameters at this energy and get the level
  // spacing and widths
  const auto parameters = this->resonanceTable()( energy );
  const auto spacing = parameters.levelSpacing();
  const Degrees degrees = this->resonanceTable().degreesOfFreedom();
  const double vl = degrees.elastic * penetrability / ratio; // ENDF D.98
  const Widths widths{ parameters.elastic() * sqrt( energy ) * vl ,
                       parameters.capture(),
                       parameters.fission(),
                       parameters.competition() };

  // calculate the fluctuation integrals
  const FluctuationIntegrals integrals =
    calculateFluctuationIntegrals( widths, degrees );

  // calculate the resulting cross sections
  result[ this->elasticID() ] +=
    factor * ( spinFactor / spacing *
    ( widths.elastic * ( widths.elastic * integrals.elastic - 2. * sin2phi ) ) );
  result[ this->captureID() ] +=
    factor * spinFactor / spacing *
    ( widths.elastic * widths.capture * integrals.capture );
  if ( widths.hasFission() ) {

    result[ this->fissionID() ] +=
      factor * spinFactor / spacing *
      ( widths.elastic * widths.fission * integrals.fission );
  }
}
