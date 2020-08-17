static
ReactionID makeReactionID( const ParticlePairID& incident,
                           const ParticlePairID& pair ) {

  if ( incident == pair ) {

    return ReactionID( incident.particle(), incident.residual(),
                       ReactionType( "elastic" ) );
  }
  return ReactionID( incident, pair );
}
