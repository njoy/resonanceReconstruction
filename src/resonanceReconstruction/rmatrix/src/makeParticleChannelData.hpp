std::vector< ParticleChannelData >
consolidateChannelData( const std::vector< ParticleChannelData >& channels ) {

  // the consolidated channels
  std::vector< ParticleChannelData > consolidated;

  // usefull lambdas
  auto isEliminated = [] ( const auto& channel ) {

    return channel.isEliminatedChannel();
  };
  auto isNotEliminated = [] ( const auto& channel ) {

    return not channel.isEliminatedChannel();
  };
  auto getEnergies = [] ( const auto& channel ) {

    return channel.energies();
  };
  auto getWidths = [] ( const auto& channel ) {

    return channel.widths();
  };
  auto getJpi = [] ( const auto& channel ) {

    return std::make_pair( channel.quantumNumbers().totalAngularMomentum(),
                           channel.quantumNumbers().parity() );
  };

  // get all data for channels that are not eliminated
  consolidated = channels | ranges::view::filter( isNotEliminated );
  auto remaining = channels | ranges::view::filter( isEliminated );

  // get the different Jpi values in the eliminated channels
  std::vector< std::pair< TotalAngularMomentum, Parity > > spins =
    channels | ranges::view::filter( isEliminated )
             | ranges::view::transform( getJpi )
             | ranges::to_vector
             | ranges::action::sort
             | ranges::action::unique;

  // go over the Jpi
  for ( const auto& Jpi : spins ) {

    auto J = Jpi.first;
    auto parity = Jpi.second;

    auto filter = [&] ( const auto& channel ) {

      return ( channel.quantumNumbers().totalAngularMomentum() == J ) and
             ( channel.quantumNumbers().parity() == parity );
    };

    auto eliminated = remaining | ranges::view::filter( filter );
    if ( ranges::distance( eliminated ) == 1 ) {

      consolidated.push_back( eliminated.front() );
    }
    else {

      auto channel = eliminated.front().channel();
      auto energies = eliminated | ranges::view::transform( getEnergies );
      auto widths = eliminated | ranges::view::transform( getWidths );
      consolidated.emplace_back( channel,
                                 std::move( energies | ranges::view::join ),
                                 std::move( widths | ranges::view::join ),
                                 true );
    }
  }

  return consolidated;
}

template < typename Formalism, typename BoundaryOption >
std::vector< ParticleChannelData >
makeParticleChannelData(
    const ParticlePair& incident,
    const std::vector< ParticlePair > pairs,
    const endf::RMatrixLimited::ParticlePairs& endfPairs,
    const endf::RMatrixLimited::SpinGroup& endfSpinGroup,
    bool reducedWidthsFlag,
    Formalism,
    BoundaryOption ) {

  std::vector< ParticleChannelData > data;

  // some usefull lambdas
  auto isIncidentChannel = [] ( const auto& channel ) {

    return channel.isIncidentChannel();
  };
  auto first = [] ( const auto& pair ) {

    return std::get< 0 >( pair );
  };
  auto second = [] ( const auto& pair ) {

    return std::get< 1 >( pair );
  };
  auto nonZero = [] ( const auto& pair ) {

    return std::get< 1 >( pair ) != 0.;
  };
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

      auto pairs = ranges::view::zip(
                     endfSpinGroup.parameters().resonanceEnergies(),
                     endfSpinGroup.parameters().resonanceParameters()
                       | ranges::view::transform( [i] ( const auto& widths )
                                                      { return widths[i]; } ) );
      auto nonzero = pairs | ranges::view::filter( nonZero );

      // only add the channel if there are resonances
      if ( ranges::distance( nonzero ) != 0 ) {

        auto energies = nonzero | ranges::view::transform( first )
                                | ranges::view::transform( toEnergy )
                                | ranges::to_vector;
        auto widths = nonzero | ranges::view::transform( second )
                              | ranges::to_vector;
        auto reduced = toReducedWidths( channels[i], energies, widths );
        data.emplace_back( channels[i], std::move( energies ),
                           std::move( reduced ), false );
      }
      else {

        // the channel has no data - only elastic channels contribute
        // if it is an elastic channel, add an empty elastic channel to take
        // into account potential scattering
        if ( std::visit( isIncidentChannel, channels[i] ) == true ) {

          data.emplace_back( channels[i],
                             std::vector< Energy >{},
                             std::vector< ReducedWidth >{},
                             false );
        }
      }
    }
    else {

      // the entire spin group is empty - only elastic channels contribute
      // if it is an elastic channel, add an empty elastic channel to take
      // into account potential scattering
      if ( std::visit( isIncidentChannel, channels[i] ) == true ) {

        data.emplace_back( channels[i],
                           std::vector< Energy >{},
                           std::vector< ReducedWidth >{},
                           false );
      }
    }
  }

  return data;
}

template < typename BoundaryOption >
std::vector< ParticleChannelData >
makeParticleChannelData(
    const ParticlePair& incident,
    const std::vector< ParticlePair > pairs,
    const endf::RMatrixLimited::ParticlePairs& endfPairs,
    const endf::RMatrixLimited::SpinGroup& endfSpinGroup,
    bool reducedWidthsFlag,
    ReichMoore,
    BoundaryOption ) {

  std::vector< ParticleChannelData > data;

  // some usefull lambdas
  auto isIncidentChannel = [] ( const auto& channel ) {

    return channel.isIncidentChannel();
  };
  auto first = [] ( const auto& pair ) {

    return std::get< 0 >( pair );
  };
  auto second = [] ( const auto& pair ) {

    return std::get< 1 >( pair );
  };
  auto nonZero = [] ( const auto& pair ) {

    return std::get< 1 >( pair ) != 0.;
  };
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

      auto pairs = ranges::view::zip(
                     endfSpinGroup.parameters().resonanceEnergies(),
                     endfSpinGroup.parameters().resonanceParameters()
                       | ranges::view::transform( [i] ( const auto& widths )
                                                      { return widths[i]; } ) );
      auto nonzero = pairs | ranges::view::filter( nonZero );

      if ( ranges::distance( nonzero ) != 0 ) {

        auto energies = nonzero | ranges::view::transform( first )
                                | ranges::view::transform( toEnergy )
                                | ranges::to_vector;
        auto widths = nonzero | ranges::view::transform( second )
                              | ranges::to_vector;
        auto reduced = toReducedWidths( channels[i], energies, widths );
        data.emplace_back( channels[i], std::move( energies ),
                           std::move( reduced ),
                           i == eliminated ? true : false );
      }
      else {

        // the channel has no data - only elastic channels contribute
        // if it is an elastic channel, add an empty elastic channel to take
        // into account potential scattering
        if ( std::visit( isIncidentChannel, channels[i] ) == true ) {

          data.emplace_back( channels[i],
                             std::vector< Energy >{},
                             std::vector< ReducedWidth >{},
                             false );
        }
      }
    }
    else {

      // the entire spin group is empty - only elastic channels contribute
      // if it is an elastic channel, add an empty elastic channel to take
      // into account potential scattering
      if ( std::visit( isIncidentChannel, channels[i] ) == true ) {

        data.emplace_back( channels[i],
                           std::vector< Energy >{},
                           std::vector< ReducedWidth >{},
                           false );
      }
    }
  }

  return data;
}
