/**
 *  @brief Evaluate the cross sections at the given energy using SLBW
 *
 *  @param[in] energy       the incident energy
 *  @param[in,out] result   a map containing the accumulated cross sections
 */
void evaluate( const Energy& energy,
               std::map< ReactionID, CrossSection >& result ) {

  // data we need: k, P, phi, rho, g_J
  const auto channel = this->incidentChannel();
  const auto incident = channel.particlePair().pairID();
  const auto incident = channel.particlePair().particle().particleID();
  const auto target = channel.particlePair().residual().particleID();
  const auto waveNumber = channel.waveNumber( energy );
  const auto qx = this->QX();
  const auto p = channel.penetrability( energy );
  const auto q = qx.value != 0. ? channel.penetrability( energy - qx ) : p;
  const auto s = channel.shiftFactor( energy );
  const auto phaseShift = channel.phaseShift( energy );
  const auto radius = channel.radii().penetrabilityRadius( energy );
  const auto spinFactor = channel.statisticalSpinFactor();
  const auto ratio = waveNumber * radius;
  const auto sinphi = std::sin( phaseShift );
  const auto sintwophi = std::sin( 2. * phaseShift );
  const auto sin2phi = sinphi * sinphi;

  // the 2 * pi2 / k2 factor
  const CrossSection factor = pi / ( waveNumber * waveNumber ) * spinFactor;

  // lambda to calculate the cross sections for each resonance
  auto calculate = [&] ( const auto& resonance ) -> Components {

    const auto total = resonance.total( p, q );
    const auto elastic = resonance.elastic( p );
    const auto capture = resonance.capture();
    const auto fission = resonance.fission();
    const auto delta = energy - resonance.energyPrime( s );
    const auto denominator = delta * delta + 0.25 * total * total;

    return { ( elastic * elastic - 2. * elastic * total * sin2phi
               + 2. * delta * total * sintwophi ) / denominator,
             capture * elastic / denominator,
             fission * elastic / denominator,
             0.0 };
  };

  // accumulate the cross section components
  const Components = ranges::accumulate(
                       this->table().resonances()
                         | ranges::view::transform( calculate ),
                       { 0., 0., 0., 0. } );

  // calculate the resulting cross sections
  ReactionID elas = ReactionID( incident, target, elementary::ReactionType( "elastic" ) );
  ReactionID capt = ReactionID( incident, target, elementary::ReactionType( "capture" ) );
  ReactionID fiss = ReactionID( incident, target, elementary::ReactionType( "fission" ) );
  result[ elas ] += factor * components.elastic;
  result[ capt ] += factor * components.capture;
  if ( components.hasFission() ) {

    result[ fiss ] += factor * components.fission;
  }
}
