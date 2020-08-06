std::vector< ChannelQuantumNumbers >
makeQuantumNumbers( const ParticlePair& incident, unsigned int lmax ) {

  std::vector< ChannelQuantumNumbers > numbers;
  auto channelspins = possibleChannelSpinValues( incident.particle().spin(),
                                                 incident.residual().spin() );
  for ( unsigned int l = 0; l < lmax; ++l ) {

    for ( unsigned int i = 0; i < channelspins.size(); ++i ) {

      auto spinvalues = possibleChannelTotalAngularMomentumValues(
                            l, channelspins[i] );
      for ( unsigned int j = 0; j < spinvalues.size(); ++j ) {

        numbers.emplace_back( l, channelspins[i], spinvalues[j],
                              std::pow( -1, l ) >= 0. ? +1 : -1 );
      }
    }
  }
  return numbers;
}
