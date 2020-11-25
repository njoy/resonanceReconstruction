/**
 *  @brief Evaluate the cross sections at the given energy using MLBW
 *
 *  @param[in] energy       the incident energy
 *  @param[in,out] result   a map containing the accumulated cross sections
 */
void evaluate( const Energy& energy,
               Map< ReactionID, CrossSection >& result ) {

  // data we need: k, P, phi, rho, g_J
  const auto channel = this->incidentChannel();
  const auto waveNumber = channel.waveNumber( energy );
  const auto spinFactor = channel.statisticalSpinFactor();

  // the pi / k2 factor
  const CrossSection factor = pi / ( waveNumber * waveNumber ) * spinFactor;

  // calculate the SLBW cross sections
  SpinGroup< SingleLevelBreitWigner >::evaluate( energy, result );

  // add the elastic cross term for MLBW
  double term = 0.;
  unsigned int nr = this->resonanceTable().resonances().size();
  auto elastic = this->elastic();
  auto total = this->total();
  auto delta = this->delta();
  auto denominator = this->denominator();
  for ( unsigned int r = 1; r < nr; ++r ) {

    for ( unsigned int rp = 0; rp < r; ++rp ) {

      term += 2. * elastic[r] * elastic[rp]
                 * ( delta[r] * delta[rp]+ 0.25 * total[r] * total[rp] )
              / denominator[r] / denominator[rp];
    }
  }

  // calculate the resonance crossterm and add it to the elastic cross section
  result[ this->elasticID() ] += factor * term;
}
