std::vector< ParticleChannelData >
makeReichMooreChannelData(
    const ENDF::resolved::ReichMooreLValue& endfLValue,
    const ParticlePair& pair, const ChannelRadii& radii,
    std::vector< ChannelQuantumNumbers >& available,
    const unsigned int& naps ) {

  // the spin groups found
  std::vector< ParticleChannelData > data;

  // useful numbers
  unsigned int l = endfLValue.orbitalMomentum();
  double apl = endfLValue.lDependentScatteringRadius();
  double awri = endfLValue.atomicWeightRatio();

  // some usefull lambdas
  auto compareNumbers = [] ( const auto& left, const auto& right ) -> bool {

    return
    ( left.orbitalAngularMomentum() == right.orbitalAngularMomentum() ) and
    ( left.spin() == right.spin() ) and
    ( left.totalAngularMomentum() == right.totalAngularMomentum() ) and
    ( left.parity() == right.parity() );
  };
  auto toEnergy = [&] ( double value ) -> Energy {

    return value * electronVolt;
  };
  auto toReducedWidth = [&] ( double value ) -> ReducedWidth {

    return value * rootElectronVolt;
  };
  auto toElasticWidth = [&] ( double value, double p ) -> ReducedWidth {

    int sign = value >= 0. ? +1 : -1;
    return toReducedWidth( sign * std::sqrt( std::abs( value ) / 2. / p ) );
  };
  auto toOtherWidth = [&] ( double value ) -> ReducedWidth {

    int sign = value >= 0. ? +1 : -1;
    return toReducedWidth( sign * std::sqrt( std::abs( value ) / 2. ) );
  };
  auto isAllZero = [&] ( auto&& range ) -> bool {

    return ranges::accumulate(
               range | ranges::view::transform(
                           [] ( const auto& value )
                              { return value == 0.; } ),
               true,
               [] ( bool left, bool right ) { return left and right; } );
  };
  auto retrieveQuantumNumber =
       [&] ( unsigned int l, double j,
             std::vector< ChannelQuantumNumbers >& available ) {

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

  // the fission and capture particle pair
  ParticlePair pFission( pair.particle(), pair.residual(),
                         ParticlePairID( "fission" ) );
  ParticlePair pCapture( pair.particle(), pair.residual(),
                         ParticlePairID( "capture" ) );

  // the channel radii to be used in the channels
  ChannelRadii useRadii = radii;
  if ( apl != 0. ) {

    useRadii = makeChannelRadii( apl, naps, awri, pair.particle().mass().value );
  }

  // the different J values referenced in the table
  std::vector< double > jvalues = endfLValue.spinValues()
                                      | ranges::to_vector
                                      | ranges::action::sort
                                      | ranges::action::unique;

  // loop over each J value
  auto resonances = endfLValue.resonances();
  for ( double J : jvalues ) {

    // filter out resonances with this J value
    auto current = resonances | ranges::view::filter(
                                     [J] ( const auto& resonance )
                                         { return resonance.spin() == J; } );

    // quantum number for this l,J pair
    ChannelQuantumNumbers nElastic = retrieveQuantumNumber( l, J, available );
    ChannelQuantumNumbers nOther( l, 0, std::abs( J ),
                                  std::pow( -1, l ) >= 0. ? +1 : -1 );

    // elastic and fission channels that will go in the spin groups
    Channel< Neutron > cElastic( pair, pair,  0.0 * electronVolt,
                                nElastic, useRadii );
    Channel< Photon > cCapture( pair, pCapture,  0.0 * electronVolt,
                                nElastic, useRadii );
    Channel< Fission > cFission1( pair, pFission, "fission1",
                                  0.0 * electronVolt,
                                  nOther, useRadii );
    Channel< Fission > cFission2( pair, pFission, "fission2",
                                  0.0 * electronVolt,
                                  nOther, useRadii );

    // get the widths
    auto wElastic = current | ranges::view::transform(
                                [] ( const auto& resonance )
                                   { return resonance.neutronWidth(); } );
    auto wCapture = current | ranges::view::transform(
                                [] ( const auto& resonance )
                                   { return resonance.gammaWidth(); } );
    auto wFission1 = current | ranges::view::transform(
                                 [] ( const auto& resonance )
                                    { return resonance.firstFissionWidth(); } );
    auto wFission2 = current | ranges::view::transform(
                                 [] ( const auto& resonance )
                                    { return resonance.secondFissionWidth(); } );

    // calculate penetrabilities
    std::vector< Energy > energies =
        current | ranges::view::transform(
                    [] ( const auto& resonance )
                       { return resonance.resonanceEnergy() * electronVolt; } );
    auto penetrabilities =
        energies | ranges::view::transform(
                     [&] ( const auto& energy )
                         { return cElastic.penetrability( energy ); } );

    // calculate reduced widths
    auto rwElastic = ranges::view::zip_with(
                         toElasticWidth, wElastic, penetrabilities );
    auto rwCapture = wCapture | ranges::view::transform( toOtherWidth );
    auto rwFission1 = wFission1 | ranges::view::transform( toOtherWidth );
    auto rwFission2 = wFission2 | ranges::view::transform( toOtherWidth );

    // see how many channels we have
    bool bFission1 = not isAllZero( wFission1 );
    bool bFission2 = not isAllZero( wFission2 );

    // add channel data
    data.emplace_back( cElastic, std::vector< Energy >( energies ), rwElastic );
    data.emplace_back( cCapture, std::vector< Energy >( energies ), rwCapture, true );
    if ( bFission1 ) {

      data.emplace_back( cFission1, std::vector< Energy >( energies ), rwFission1 );
    }
    if ( bFission2 ) {

      if ( bFission1 ) {

        data.emplace_back( cFission2, std::vector< Energy >( energies ), rwFission2 );
      }
      else {

        data.emplace_back( cFission1, std::vector< Energy >( energies ), rwFission2 );
      }
    }
  }

  return data;
}
