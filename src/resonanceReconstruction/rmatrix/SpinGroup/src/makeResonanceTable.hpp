static
ResonanceTable
makeResonanceTable( const std::vector< ParticleChannelData >& channels ) {

  auto getID = [] ( const auto& data )
                  { return data.channelID(); };
  auto isEliminated = [] ( const auto& data )
                         { return data.isEliminatedChannel(); };
  auto isTrue = [] ( const auto& value )
                   { return value == true; };

  // there may be at most one eleminated channel
  auto eliminated = channels | ranges::views::transform( isEliminated );
  auto indexEliminated = ranges::cpp20::distance(
                           ranges::cpp20::begin( eliminated ),
                           std::find_if( ranges::cpp20::begin( eliminated ),
                                         ranges::cpp20::end( eliminated ),
                                        isTrue ) );
  auto numberEliminated = ranges::cpp20::distance(
                              ranges::cpp20::views::filter( eliminated, isTrue ) );
  if ( numberEliminated > 1 ) {

    Log::error( "More than 1 eliminated channel was found." );
    throw std::exception();
  }

  // get the channel IDs (remove the eliminated channel if required)
  std::vector< ChannelID > ids =
      ranges::to< std::vector< ChannelID > >(
          channels | ranges::views::transform( getID ) );
  auto numberChannels = ranges::cpp20::distance( ids );
  if ( numberEliminated ) {

    ids.erase( ids.begin() + indexEliminated );
  }

  // get the widths by energy
  std::map< Energy, std::vector< ReducedWidth > > map;
  auto addWidth = [&] ( const auto& energy, const auto& width, unsigned int index ) {

    map.try_emplace( energy, numberChannels ).first->second[ index ] = width;
  };
  auto addWidthsToMap = [&] ( auto&& tuple ) {

    auto data = std::get<0>( tuple );
    auto index = std::get<1>( tuple );

    auto energies = data.energies();
    auto widhts = data.widths();
    auto indices = ranges::views::repeat_n( index,
                                            ranges::cpp20::distance( energies ) );

    ranges::cpp20::for_each( ranges::views::zip( energies, widhts, indices ),
                             [&] ( auto&& tuple )
                                 { addWidth( std::get<0>( tuple ),
                                             std::get<1>( tuple ),
                                             std::get<2>( tuple ) ); } );
  };
  ranges::cpp20::for_each( ranges::views::zip(
                               channels,
                               ranges::cpp20::views::iota( 0, numberChannels ) ),
                           addWidthsToMap );

  // create the resonances
  auto toResonance = [&] ( const auto& pair ) {

    auto energy = pair.first;
    auto widths = pair.second;
    ReducedWidth eWidth = numberEliminated ? widths[ indexEliminated ] :
                                             0. * rootElectronVolt;

    if ( numberEliminated ) {

      widths.erase( widths.begin() + indexEliminated );
    }

    return Resonance( energy, std::move( widths ), eWidth );
  };
  std::vector< Resonance > resonances =
      ranges::to< std::vector< Resonance > >(
          map | ranges::views::transform( toResonance ) );

  return ResonanceTable( std::move( ids ), std::move( resonances ) );
}
