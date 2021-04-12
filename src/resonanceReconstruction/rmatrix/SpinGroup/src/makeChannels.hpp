static
std::vector< ParticleChannel >
makeChannels( const std::vector< ParticleChannelData >& channels ) {

  std::vector< ParticleChannel > result =
      ranges::to< std::vector< ParticleChannel > >(
          channels | ranges::views::transform( [] ( const auto& data )
                                                  { return data.channel(); } ) );

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
  if ( indexEliminated != ranges::cpp20::distance( channels ) ) {

    result.erase( result.begin() + indexEliminated );
  }

  return result;
}
