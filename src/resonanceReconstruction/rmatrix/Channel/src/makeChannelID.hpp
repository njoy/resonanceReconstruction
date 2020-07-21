static
ChannelID makeChannelID( const ParticlePairID& id,
                         const ChannelQuantumNumbers& numbers ) {

  return id.symbol() + numbers.toString();
}
