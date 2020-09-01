template < typename Formalism, typename BoundaryOption >
std::vector< ParticleChannelData >
makeParticleChannelData(
    const ParticlePair& incident,
    const std::vector< ParticlePair > pairs,
    const ENDF::resolved::RMatrixLimited::ParticlePairs& endfPairs,
    const ENDF::resolved::RMatrixLimited::SpinGroup& endfSpinGroup,
    bool reducedWidthsFlag,
    Formalism,
    BoundaryOption ) {

  std::vector< ParticleChannelData > data;

  // some usefull lambdas
  auto toEnergy = [&] ( double value ) -> Energy {

    return value * electronVolt;
  };
  auto toReducedWidth = [&] ( double value ) -> ReducedWidth {

    return value * rootElectronVolt;
  };
  auto calculateReducedWidth = [&] ( const Energy& energy,
                                     double width,
                                     const ParticleChannel& channel ) {

    double P = std::visit( [&] ( const auto& channel )
                               { return channel.penetrability( energy ); },
                             channel );
    return toReducedWidth( ( width < 0. ? -1. : +1. ) *
                           std::sqrt( std::abs( width ) / 2. / P ) );
  };
  auto toReducedWidths = [&] ( const ParticleChannel& channel,
                               const auto& energies,
                               const auto& widths ) -> std::vector< ReducedWidth > {

    if ( reducedWidthsFlag ) {

      return widths | ranges::view::transform( toReducedWidth );
    }

    return ranges::view::zip_with(
               calculateReducedWidth,
               energies,
               widths,
               ranges::view::repeat_n( channel, energies.size() ) );
  };

  std::vector< ParticleChannel > channels =
      makeParticleChannels( incident, pairs, endfPairs,
                            endfSpinGroup.channels() );

  // go over the data and create the particle channel data
  for ( unsigned int i = 0; i < channels.size(); ++i ) {

    if ( endfSpinGroup.parameters().numberResonances() != 0 ) {

      auto energies = endfSpinGroup.parameters().resonanceEnergies()
                        | ranges::view::transform( toEnergy );
      auto widths = endfSpinGroup.parameters().resonanceParameters()
                      | ranges::view::transform( [i] ( const auto& widths )
                                                     { return widths[i]; } );
      auto reduced = toReducedWidths( channels[i], energies, widths );
      data.emplace_back( channels[i], std::move( energies ),
                         std::move( reduced ), false );
    }
    else {

      data.emplace_back( channels[i],
                         std::vector< Energy >{},
                         std::vector< ReducedWidth >{},
                         false );
    }
  }

  return data;
}

template < typename BoundaryOption >
std::vector< ParticleChannelData >
makeParticleChannelData(
    const ParticlePair& incident,
    const std::vector< ParticlePair > pairs,
    const ENDF::resolved::RMatrixLimited::ParticlePairs& endfPairs,
    const ENDF::resolved::RMatrixLimited::SpinGroup& endfSpinGroup,
    bool reducedWidthsFlag,
    ReichMoore,
    BoundaryOption ) {

  std::vector< ParticleChannelData > data;

  // some usefull lambdas
  auto toEnergy = [&] ( double value ) -> Energy {

    return value * electronVolt;
  };
  auto toReducedWidth = [&] ( double value ) -> ReducedWidth {

    return value * rootElectronVolt;
  };
  auto calculateReducedWidth = [&] ( const Energy& energy,
                                     double width,
                                     const ParticleChannel& channel ) {

    double P = std::visit( [&] ( const auto& channel )
                               { return channel.penetrability( energy ); },
                             channel );
    return toReducedWidth( ( width < 0. ? -1. : +1. ) *
                           std::sqrt( std::abs( width ) / 2. / P ) );
  };
  auto toReducedWidths = [&] ( const ParticleChannel& channel,
                               const auto& energies,
                               const auto& widths ) -> std::vector< ReducedWidth > {

    if ( reducedWidthsFlag ) {

      return widths | ranges::view::transform( toReducedWidth );
    }

    return ranges::view::zip_with(
               calculateReducedWidth,
               energies,
               widths,
               ranges::view::repeat_n( channel, energies.size() ) );
  };

  std::vector< ParticleChannel > channels =
      makeParticleChannels( incident, pairs, endfPairs,
                            endfSpinGroup.channels() );
  unsigned int eliminated = [&] {

    auto getParticlePair = [] ( const ParticleChannel& channel ) {

      return std::visit( [] ( const auto& channel )
                            { return channel.particlePair(); },
                         channel );
    };

    auto eliminatedPair = pairs[ rmatrix::eliminated( endfPairs ) ];
    auto pairs = channels | ranges::view::transform( getParticlePair );
    auto found = std::find_if(
                     ranges::begin( pairs ),
                     ranges::end( pairs ),
                     [&] ( const auto& pair )
                         { return eliminatedPair.pairID() == pair.pairID(); } );
    return std::distance( ranges::begin( pairs ), found );
  }();

  // go over the data and create the particle channel data
  for ( unsigned int i = 0; i < channels.size(); ++i ) {

    if ( endfSpinGroup.parameters().numberResonances() != 0 ) {

      auto energies = endfSpinGroup.parameters().resonanceEnergies()
                        | ranges::view::transform( toEnergy );
      auto widths = endfSpinGroup.parameters().resonanceParameters()
                      | ranges::view::transform( [i] ( const auto& widths )
                                                     { return widths[i]; } );
      auto reduced = toReducedWidths( channels[i], energies, widths );
      data.emplace_back( channels[i], std::move( energies ),
                         std::move( reduced ), i == eliminated ? true : false );
    }
    else {

      data.emplace_back( channels[i],
                         std::vector< Energy >{},
                         std::vector< ReducedWidth >{},
                         i == eliminated ? true : false );
    }
  }

  return data;
}
