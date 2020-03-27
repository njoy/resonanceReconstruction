inline std::vector< ParticleChannel >
makeParticleChannels(
    const std::vector< ParticlePair >& pairs,
    const ENDF::resolved::RMatrixLimited::ParticlePairs& endfPairs,
    const ENDF::resolved::RMatrixLimited::ResonanceChannels& endfChannels ) {

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
  auto makeParticleChannel = [&] ( const ParticlePair& pair,
                                   const ChannelQuantumNumbers& numbers,
                                   const ChannelRadii& radii,
                                   const BoundaryCondition& boundary,
                                   int mt ) -> ParticleChannel {

    switch ( mt ) {

      case 102 : return Channel< Photon >{ pair, numbers, radii, boundary };
      case  18 : return Channel< Fission >{ pair, numbers, radii, boundary };
      case  19 : return Channel< Fission >{ pair, numbers, radii, boundary };
      case  20 : return Channel< Fission >{ pair, numbers, radii, boundary };
      case  21 : return Channel< Fission >{ pair, numbers, radii, boundary };
      case  38 : return Channel< Fission >{ pair, numbers, radii, boundary };
      default : {

        if ( pair.particle().charge() > 0. * coulomb ) {

          return Channel< ChargedParticle >{ pair, numbers, radii, boundary };
        }
        else {

          return Channel< Neutron >{ pair, numbers, radii, boundary };
        }
      }
    }
  };

  // do some range magic
  auto channelPairs = endfChannels.particlePairNumbers()
                          | ranges::view::transform( getParticlePair );
  auto channelNumbers = ranges::view::zip_with(
                            makeChannelQuantumNumbers,
                            endfChannels.orbitalMomentumValues(),
                            endfChannels.channelSpinValues() );
  auto channelRadii = ranges::view::zip_with(
                          makeChannelRadii,
                          endfChannels.trueChannelRadii(),
                          endfChannels.effectiveChannelRadii() );

  return ranges::view::zip_with(
             makeParticleChannel,
             channelPairs,
             channelNumbers,
             channelRadii,
             endfChannels.boundaryConditionValues(),
             endfPairs.MT() );
}
