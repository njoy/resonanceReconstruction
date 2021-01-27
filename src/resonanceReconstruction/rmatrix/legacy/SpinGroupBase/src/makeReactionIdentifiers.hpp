static
std::vector< ReactionID >
makeReactionIdentifiers( const Channel< Neutron >& channel,
                         const ResonanceTableType& table ) {

  std::vector< ReactionID > reactions;

  // add elastic and capture
  auto incident = channel.particlePair().particle().particleID();
  auto target = channel.particlePair().residual().particleID();
  reactions.push_back( ReactionID{ incident, target, ReactionType( "elastic" ) } );
  reactions.push_back( ReactionID{ incident, target, ReactionType( "capture" ) } );

  // lambdas
  auto hasFission = [] ( const auto& resonance )
                       { return resonance.hasFission(); };

  // add fission if it exists
  if ( ranges::count_if( table.resonances(), hasFission ) > 1 ) {

    reactions.push_back( ReactionID{ incident, target, ReactionType( "fission" ) } );
  }

  return reactions;
}
