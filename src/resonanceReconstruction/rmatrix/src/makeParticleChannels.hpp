std::vector< ParticleChannel >
makeParticleChannels(
    const ParticlePair& incident,
    const std::vector< ParticlePair >& pairs,
    const endf::RMatrixLimited::ParticlePairs& endfPairs,
    const endf::RMatrixLimited::ResonanceChannels& endfChannels ) {

  // ENDF information we need
  double spin = endfChannels.spin();
  int parity = endfChannels.parity();
  double j = std::abs( spin );
  Parity pi = spin == 0.0 ? ( parity >= 0 ? Parity( +1 ) : Parity( -1 ) )
                          : ( spin > 0. ? Parity( +1 ) : Parity( -1 ) );

  // a few useful lambdas
  auto getParticlePair = [&] ( unsigned int i ) { return pairs[ i - 1 ]; };
  auto makeChannelQuantumNumbers = [&] ( unsigned int l, double s ) {

    return ChannelQuantumNumbers{ l, s, j, pi };
  };
  auto makeChannelRadii = [&] ( double apt, double ape ) {

    return ChannelRadii{ apt * rootBarn, ape * rootBarn };
  };
  auto makeParticleChannel = [&] ( const ParticlePair& incident,
                                   const ParticlePair& pair,
                                   double Q,
                                   const ChannelQuantumNumbers& numbers,
                                   const ChannelRadii& radii,
                                   const BoundaryCondition& boundary,
                                   int mt ) -> ParticleChannel {

    switch ( mt ) {

      case 102 : return Channel< Photon >{ incident, pair, Q * electronVolt,
                                           numbers, radii, boundary };
      case  18 : return Channel< Fission >{ incident, pair, Q * electronVolt,
                                            numbers, radii, boundary };
      case  19 : return Channel< Fission >{ incident, pair, Q * electronVolt,
                                            numbers, radii, boundary };
      case  20 : return Channel< Fission >{ incident, pair, Q * electronVolt,
                                            numbers, radii, boundary };
      case  21 : return Channel< Fission >{ incident, pair, Q * electronVolt,
                                            numbers, radii, boundary };
      case  38 : return Channel< Fission >{ incident, pair, Q * electronVolt,
                                            numbers, radii, boundary };
      default : {

        if ( pair.particle().charge() > 0. * coulomb ) {

          return Channel< ChargedParticle >{ incident, pair, Q * electronVolt,
                                             numbers, radii, boundary };
        }
        else {

          return Channel< Neutron >{ incident, pair, Q * electronVolt,
                                     numbers, radii, boundary };
        }
      }
    }
  };

  // do some range magic
  auto qPairs = endfPairs.Q();
  auto mtPairs = endfPairs.MT();
  auto incidentPairs = ranges::views::repeat_n( incident,
                                                endfChannels.numberChannels() );
  auto channelPairs = endfChannels.particlePairNumbers()
                          | ranges::cpp20::views::transform( getParticlePair );
  auto qValues = endfChannels.particlePairNumbers()
                     | ranges::cpp20::views::transform(
                           [&] ( unsigned int i )
                               { return qPairs[ i - 1 ]; } );
  auto channelNumbers = ranges::views::zip_with(
                            makeChannelQuantumNumbers,
                            endfChannels.orbitalMomentumValues(),
                            endfChannels.channelSpinValues() );
  auto channelRadii = ranges::views::zip_with(
                          makeChannelRadii,
                          endfChannels.trueChannelRadii(),
                          endfChannels.effectiveChannelRadii() );
  auto mtNumbers = endfChannels.particlePairNumbers()
                       | ranges::views::transform(
                             [&] ( unsigned int i )
                                 { return mtPairs[ i - 1 ]; } );

  return ranges::to< std::vector< ParticleChannel > >(
           ranges::views::zip_with(
             makeParticleChannel,
             incidentPairs,
             channelPairs,
             qValues,
             channelNumbers,
             channelRadii,
             endfChannels.boundaryConditionValues(),
             mtNumbers ) );
}
