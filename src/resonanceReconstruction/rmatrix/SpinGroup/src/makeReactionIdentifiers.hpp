static
std::vector< ReactionID >
makeReactionIdentifiers( const std::vector< ParticleChannel >& channels,
                         unsigned int incident ) {

  const unsigned size = channels.size();
  if ( incident > size - 1 ) {

    Log::error( "Erroneous incident channel index." );
    Log::info( "Expected an index less than or equal to {}, found {}",
               size - 1, incident );
    throw std::exception();
  }

  const auto pairs =
    channels
      | ranges::view::transform(
            [&] ( const auto& channel )
                { return std::visit(
                      [&] ( const auto& channel )
                          { return channel.particlePair().pairID(); },
                      channel ); } );
  const auto in = pairs[ incident ];
  return ranges::view::concat(
             pairs | ranges::view::transform(
                         [&] ( const auto& pair )
                             { return ReactionID( in + "->" + pair ); } ),
             ranges::view::single( in + "->capture" ) );
}
