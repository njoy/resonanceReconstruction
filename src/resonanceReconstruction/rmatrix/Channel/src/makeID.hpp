static
ChannelID makeID( const ParticlePairID& id,
                  const ChannelQuantumNumbers& numbers ) {
  return id + numbers.toString();
}

