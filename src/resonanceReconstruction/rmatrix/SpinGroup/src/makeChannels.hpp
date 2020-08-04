static
std::vector< ParticleChannel >
makeChannels( const std::vector< ParticleChannelData >& channels ) {

  std::vector< ParticleChannel > result =
      channels | ranges::view::transform( [] ( const auto& data )
                                             { return data.channel(); } );

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
  if ( indexEliminated != ranges::distance( channels ) ) {

    result.erase( result.begin() + indexEliminated );
  }

  return result;
}
