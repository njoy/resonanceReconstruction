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

ChannelQuantumNumbers
retrieveQuantumNumber( unsigned int l, double j,
                       std::vector< ChannelQuantumNumbers >& available ) {

  // useful lambdas
  auto compareNumbers = [] ( const auto& left, const auto& right ) -> bool {

    return
    ( left.orbitalAngularMomentum() == right.orbitalAngularMomentum() ) and
    ( left.spin() == right.spin() ) and
    ( left.totalAngularMomentum() == right.totalAngularMomentum() ) and
    ( left.parity() == right.parity() );
  };

  // retrieve the available numbers with this l,J
  auto filter =
       [l,j] ( const auto& number )
             { return ( number.orbitalAngularMomentum() == l ) and
                      ( number.totalAngularMomentum() == std::abs( j ) ); };
  auto filtered = available | ranges::view::filter( filter );

  auto found = ranges::distance( filtered );
  if ( found > 0 ) {

    // retrieve the one we need and erase it from the available numbers
    ChannelQuantumNumbers numbers( 0, 0., 0., +1 );
    if ( found == 1 ) {

      numbers = filtered.front();
    }
    else {

      numbers =  filtered.front().spin() < filtered.back().spin() ?
                 filtered.front() : filtered.back();
    }
    available.erase( std::find_if(
                         available.begin(), available.end(),
                         [&] ( const auto& value )
                             { return compareNumbers( value, numbers ); } ) );
    return numbers;
  }
  else {

    throw std::runtime_error( "None of the expected spin groups has l="
                              + std::to_string( l ) + " and J="
                              + std::to_string( j ) );
  }
};
