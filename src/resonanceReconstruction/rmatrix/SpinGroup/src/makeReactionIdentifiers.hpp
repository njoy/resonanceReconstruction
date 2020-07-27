static
std::vector< ReactionID >
makeReactionIdentifiers( const std::vector< ParticleChannel >& channels,
                         ReichMoore ) {

  auto reactionID = [] ( const auto& channel )
                       { return channel.reactionID(); };

  const auto in = std::visit(
                    [] ( const auto& channel )
                       { return channel.incidentParticlePair().pairID(); },
                    channels.front() );
  return ranges::view::concat(
             channels
               | ranges::view::transform(
                     [&] ( const auto& channel )
                         { return std::visit( reactionID, channel ); } ),
             ranges::view::single( ReactionID( in, ParticlePairID( "capture" ) ) ) );
}
