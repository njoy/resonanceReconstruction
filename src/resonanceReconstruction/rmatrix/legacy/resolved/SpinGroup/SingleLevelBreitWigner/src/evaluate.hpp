/**
 *  @brief Evaluate the cross sections at the given energy using SLBW
 *
 *  @param[in] energy       the incident energy
 *  @param[in,out] result   a map containing the accumulated cross sections
 */
void evaluate( const Energy& energy,
               Map< ReactionID, CrossSection >& result ) {

  // data we need
  decltype(auto) channel = this->incidentChannel();
  const auto waveNumber = channel.waveNumber( energy );
  const auto phaseShift = channel.phaseShift( energy );
  const auto spinFactor = channel.statisticalSpinFactor();
  const auto sinphi = std::sin( phaseShift );
  const auto sintwophi = std::sin( 2. * phaseShift );
  const auto sin2phi = sinphi * sinphi;

  // the pi / k2 factor
  const CrossSection factor = pi / ( waveNumber * waveNumber ) * spinFactor;

  // precompute values for SLBW and MLBW
  this->precompute( energy );

  // lambda to calculate the cross sections for each resonance
  auto calculate = [&] ( const auto& elastic, const auto& capture,
                         const auto& fission, const auto& total,
                         const auto& delta, const auto& denominator )
                       -> Data< double > {

    return { ( elastic * ( elastic - 2. * total * sin2phi
                           + 2. * delta * sintwophi ) ) / denominator,
             capture * elastic / denominator,
             fission * elastic / denominator,
             0.0 };
  };

  // accumulate the cross section components
  const Data< double > components =
    ranges::accumulate( ranges::views::zip_with(
                            calculate,
                            this->elastic(), this->capture(), this->fission(),
                            this->total(), this->delta(), this->denominator() ),
                        Data< double >{ 0., 0., 0., 0. } );

  // calculate the resulting cross sections
  result[ this->elasticID() ] += factor * components.elastic;
  result[ this->captureID() ] += factor * components.capture;
  if ( this->hasFission() ) {

    result[ this->fissionID() ] += factor * components.fission;
  }
}
