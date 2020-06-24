static
ChannelID makeChannelID( const ParticlePairID& id,
                         const ChannelQuantumNumbers& numbers ) {

  return id + numbers.toString();
}
