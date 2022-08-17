static
std::vector< ReactionID >
makeReactionIdentifiers( const std::vector< ParticleChannel >& channels,
                         ReichMoore ) {

  std::vector< ReactionID > identifiers;
  for ( const auto& channel : channels ) {

    identifiers.emplace_back(
        std::visit( [] ( const auto& channel )
                       { return channel.reactionID(); },
                    channel ) );
  }

  const auto in = std::visit(
                    [] ( const auto& channel )
                       { return channel.incidentParticlePair().pairID(); },
                    channels.front() );

  identifiers.emplace_back( ReactionID( in, ParticlePairID( "capture" ) ) );

  return identifiers;
}

static
std::vector< ReactionID >
makeReactionIdentifiers( const std::vector< ParticleChannel >& channels,
                         GeneralRMatrix ) {

  std::vector< ReactionID > identifiers;
  for ( const auto& channel : channels ) {

    identifiers.emplace_back(
        std::visit( [] ( const auto& channel )
                       { return channel.reactionID(); },
                    channel ) );
  }
  return identifiers;
}
