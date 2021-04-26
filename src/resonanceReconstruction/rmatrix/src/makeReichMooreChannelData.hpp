std::vector< ParticleChannelData >
makeReichMooreChannelData(
    const endf::ReichMooreLValue& endfLValue,
    const AtomicMass& neutronMass,
    const ParticleID& incident,
    const ParticleID& target,
    std::vector< ChannelQuantumNumbers >& available,
    double spin,
    double ap,
    const std::optional< ChannelRadiusTable >& nro,
    const unsigned int& naps ) {

  // the spin groups found
  std::vector< ParticleChannelData > data;

  // useful numbers
  unsigned int l = endfLValue.orbitalMomentum();
  double apl = endfLValue.lDependentScatteringRadius();
  double awri = endfLValue.atomicWeightRatio();

  // some usefull lambdas
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

  // the elastic, fission and capture particle pair
  ParticlePair in( Particle( incident, neutronMass,
                             0.0 * coulombs, 0.5, +1),
                   Particle( target,
                             awri * neutronMass,
                             0.0 * coulombs, spin, +1) );
  ParticlePair pFission( in.particle(), in.residual(),
                         ParticlePairID( "fission" ) );
  ParticlePair pCapture( in.particle(), in.residual(),
                         ParticlePairID( "capture" ) );

  // the channel radii to be used in the channels
  ChannelRadii radii = makeChannelRadii( apl != 0. ? apl : ap, nro, naps,
                                         awri, neutronMass.value );

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
                                  std::pow( -1, l ) >= 0.
                                      ? static_cast< Parity >( +1 )
                                      : static_cast< Parity >( -1 ) );

    // elastic and fission channels that will go in the spin groups
    Channel< Neutron > cElastic( in, in,  0.0 * electronVolt,
                                 nElastic, radii );
    Channel< Photon > cCapture( in, pCapture,  0.0 * electronVolt,
                                nOther, radii );
    Channel< Fission > cFission1( in, pFission, "fission1",
                                  0.0 * electronVolt,
                                  nOther, radii );
    Channel< Fission > cFission2( in, pFission, "fission2",
                                  0.0 * electronVolt,
                                  nOther, radii );

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
    bool bElastic = not isAllZero( wElastic );
    bool bCapture = not isAllZero( wCapture );
    bool bFission1 = not isAllZero( wFission1 );
    bool bFission2 = not isAllZero( wFission2 );

    // add channel data
    if ( bElastic ) {

      data.emplace_back( cElastic, std::vector< Energy >( energies ), rwElastic );
    }
    else {

      data.emplace_back( cElastic, std::vector< Energy >(), std::vector< ReducedWidth >() );
    }
    if ( bCapture ) {

      data.emplace_back( cCapture, std::vector< Energy >( energies ), rwCapture, true );
    }
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
