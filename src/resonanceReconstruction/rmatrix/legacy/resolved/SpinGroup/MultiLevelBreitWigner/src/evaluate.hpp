/**
 *  @brief Evaluate the cross sections at the given energy using MLBW
 *
 *  @param[in] energy       the incident energy
 *  @param[in,out] result   a map containing the accumulated cross sections
 */
void evaluate( const Energy& energy,
               std::map< ReactionID, CrossSection >& result ) {

  // data we need: k, P, phi, rho, g_J
  const auto channel = this->incidentChannel();
  const auto waveNumber = channel.waveNumber( energy );
  const auto qx = this->QX();
  const auto p = channel.penetrability( energy );
  const auto q = qx.value != 0. ? channel.penetrability( energy - qx ) : p;
  const auto s = channel.shiftFactor( energy );
  const auto spinFactor = channel.statisticalSpinFactor();

  // the pi / k2 factor
  const CrossSection factor = pi / ( waveNumber * waveNumber ) * spinFactor;

  // resonance information
  auto resonances = this->resonanceTable().resonances();
  auto total = resonances.front().total( p, q );
  auto elastic = resonances.front().elastic( p );
  auto delta = energy - resonances.front().energyPrime( s );
  auto denominator = delta * delta + 0.25 * total * total;

  double term = 0.;
  for ( unsigned int i = 1; i < resonances.size(); ++i ) {

    auto resonancep = resonances[i];
    auto totalp = resonancep.total( p, q );
    auto elasticp = resonancep.elastic( p );
    auto deltap = energy - resonancep.energyPrime( s );
    auto denominatorp = deltap * deltap + 0.25 * totalp * totalp;

    term += 2. * elastic * elasticp * ( delta * deltap + 0.25 * total * total )
            / denominator / denominatorp;

    total = totalp;
    elastic = elasticp;
    delta = deltap;
    denominator = denominatorp;
  }

  // calculate the SLBW cross sections
  SpinGroup< SingleLevelBreitWigner >::evaluate( energy, result );

  // calculate the resonance crossterm and add it to the elastic cross section
  result[ this->elastic() ] += factor * term;
}
