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
  auto eliminated = channels | ranges::view::transform( isEliminated );
  auto indexEliminated = ranges::distance(
                           ranges::begin( eliminated ),
                           std::find_if( ranges::begin( eliminated ),
                                         ranges::end( eliminated ),
                                        isTrue ) );
  auto numberEliminated = ranges::distance(
                              ranges::view::filter( eliminated, isTrue ) );
  if ( numberEliminated > 1 ) {

    Log::error( "More than 1 eliminated channel was found." );
    throw std::exception();
  }

  // get the channel IDs (remove the eliminated channel if required)
  std::vector< ChannelID > ids = channels | ranges::view::transform( getID );
  auto numberChannels = ranges::distance( ids );
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
    auto indices = ranges::view::repeat_n( index, ranges::distance( energies ) );

    ranges::for_each( ranges::view::zip( energies, widhts, indices ),
                      [&] ( auto&& tuple )
                          { addWidth( std::get<0>( tuple ),
                                      std::get<1>( tuple ),
                                      std::get<2>( tuple ) ); } );
  };
  ranges::for_each( ranges::view::zip(
                        channels,
                        ranges::view::iota( 0, numberChannels ) ),
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
  std::vector< Resonance > resonances = map | ranges::view::transform( toResonance );

  return ResonanceTable( std::move( ids ), std::move( resonances ) );
}
