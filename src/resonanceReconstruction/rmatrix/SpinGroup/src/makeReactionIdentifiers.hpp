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
  return ranges::to< std::vector< ReactionID > >(
             ranges::views::concat(
                 channels
                   | ranges::views::transform(
                         [&] ( const auto& channel )
                             { return std::visit( reactionID, channel ); } ),
                 ranges::views::single(
                     ReactionID( in, ParticlePairID( "capture" ) ) ) ) );
}
