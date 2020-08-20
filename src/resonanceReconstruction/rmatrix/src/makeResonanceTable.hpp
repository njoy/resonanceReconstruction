ResonanceTable
makeResonanceTable(
    const std::vector< ParticleChannel > channels,
    const ENDF::resolved::RMatrixLimited::ResonanceParameters& endfParameters,
    bool reducedWidthsFlag,
    int eliminated = -1 ) {

  // usefull numbers
  unsigned int number = channels.size();

  // some usefull lambdas
  auto getID = [&] ( const ParticleChannel& channel ) {

    return std::visit( [] ( const auto& channel )
                          { return channel.channelID(); },
                       channel );
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
  auto toReducedWidths = [&] ( const Energy& energy,
                               const auto& widths ) -> std::vector< ReducedWidth > {

    if ( reducedWidthsFlag ) {

      return widths | ranges::view::transform( toReducedWidth );
    }

    return ranges::view::zip_with(
               calculateReducedWidth,
               ranges::view::repeat_n( energy, number ),
               widths,
               channels );
  };
  auto toResonance = [&] ( const Energy& energy,
                           const auto& widths ) {

    std::vector< ReducedWidth > reduced =
        toReducedWidths( energy, widths | ranges::view::take_exactly( number ) );
    if ( eliminated >= 0 ) {

      ReducedWidth eliminatedWidth = reduced[ eliminated ];
      reduced.erase( reduced.begin() + eliminated );
      return Resonance( energy, std::move( reduced ), eliminatedWidth );
    }
    return Resonance( energy, std::move( reduced ) );
  };

  // do some ranges magic
  std::vector< ChannelID > id = channels | ranges::view::transform( getID );
  if ( eliminated >= 0 ) {

    id.erase( id.begin() + eliminated );
  }
  if ( endfParameters.numberResonances() ) {

    auto resonances =
        ranges::view::zip_with( toResonance,
                                endfParameters.resonanceEnergies()
                                  | ranges::view::transform( toEnergy ),
                                endfParameters.resonanceParameters() );

    return ResonanceTable( std::move( id ), std::move( resonances ) );
  }
  else {

    return ResonanceTable( std::move( id ), {} );
  }
}
