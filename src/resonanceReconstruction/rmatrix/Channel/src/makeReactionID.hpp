static
ReactionID makeReactionID( const ParticlePairID& in,
                           const ParticlePairID& out ) {

  return in + "->" + out;
}
