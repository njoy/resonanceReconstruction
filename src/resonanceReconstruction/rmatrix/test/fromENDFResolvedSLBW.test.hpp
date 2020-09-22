std::string Cf252();

SCENARIO( "fromENDF - LRF1" ) {

  GIVEN( "valid ENDF data for Cf252" ) {

    std::string string = Cf252();
    auto begin = string.begin();
    auto end = string.end();
    long lineNumber = 1;

    njoy::ENDFtk::HeadRecord head( begin, end, lineNumber );
    njoy::ENDFtk::section::Type< 2, 151 > endf( head, begin, end, lineNumber, 9861 );
    ResonanceRange endfResonanceRange = endf.isotopes().front().resonanceRanges().front();

    auto resonances = fromENDF( endfResonanceRange, neutronMass, elementaryCharge, ParticleID( "n" ), ParticleID( "Cf252" ) );

    THEN( "the appropriate CompoundSystem is returned" ) {

      CHECK( true == resonances.isResolved() );
      CHECK( false == resonances.isUnresolved() );
      CHECK( 1e-5 == Approx( resonances.lowerEnergy().value ) );
      CHECK( 366.5 == Approx( resonances.upperEnergy().value ) );

      auto compoundsystem = std::get< legacy::resolved::CompoundSystem< SingleLevelBreitWigner > >( resonances.compoundSystem() );

      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      // content verification
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      // spin groups
      auto spingroups = compoundsystem.spinGroups();
      CHECK( 1 == spingroups.size() );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0 -  empty spin group
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup0 = spingroups[0];
/*
      // channels
      auto channels0 = spingroup0.channels();

      CHECK( 2 == channels0.size() ); // 2 normal channel + 1 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel00 = std::get< Channel< Neutron > >( channels0[0] );
      CHECK( "n,Cf252->n,Cf252" == channel00.reactionID().symbol() );

      // incident particle pair
      const auto incident00 = channel00.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident00.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident00.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident00.particle().spin() ) );
      CHECK( +1 == incident00.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( incident00.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( incident00.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident00.residual().spin() ) );
      CHECK( +1 == incident00.residual().parity() );
      CHECK( "n,Cf252" == incident00.pairID().symbol() );

      // particle pair
      const auto pair00 = channel00.particlePair();
      CHECK( 1.008664 == Approx( pair00.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair00.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair00.particle().spin() ) );
      CHECK( +1 == pair00.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( pair00.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( pair00.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair00.residual().spin() ) );
      CHECK( +1 == pair00.residual().parity() );
      CHECK( "n,Cf252" == pair00.pairID().symbol() );

      // quantum numbers
      const auto numbers00 = channel00.quantumNumbers();
      CHECK( 1 == numbers00.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers00.spin() );
      CHECK( 0.0 == numbers00.totalAngularMomentum() );
      CHECK( -1 == numbers00.parity() );
      CHECK( "{1,1,0-}" == numbers00.toString() );

      // radii
      const auto radii00 = channel00.radii();
      CHECK( .9410000 == Approx( radii00.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii00.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii00.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel00.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel00.Q().value ) );

      // resonance table
      auto table1 = spingroup1.resonanceTable();

      CHECK( 3 == table1.numberChannels() ); // 1 normal channel + 1 eliminated
      CHECK( 231 == table1.numberResonances() );

      auto energies1 = table1.energies();
      CHECK( -6.908700 == Approx( energies1.front().value ) );
      CHECK( 2.511660e+3 == Approx( energies1.back().value ) );

      auto resonances1 = table1.resonances();
      CHECK( -6.908700 == Approx( resonances1.front().energy().value ) );
      CHECK( 2.511660e+3 == Approx( resonances1.back().energy().value ) );
      CHECK( 3 == resonances1.front().widths().size() );
      CHECK( 3 == resonances1.back().widths().size() );
      CHECK( std::sqrt( 1.805854e-2 / 2. / channel10.penetrability( -6.908700 * electronVolt ) )
             == Approx( resonances1.front().widths()[0].value ) );
      CHECK( std::sqrt( 3.123558e-1 / 2. / channel10.penetrability( 2.511660e+3 * electronVolt ) )
             == Approx( resonances1.back().widths()[0].value ) );
      CHECK( -std::sqrt( 9.235517e-1 / 2. )
             == Approx( resonances1.front().widths()[1].value ) );
      CHECK( std::sqrt( 1.002582 / 2. )
             == Approx( resonances1.back().widths()[1].value ) );
      CHECK( std::sqrt( 3.464891e-1 / 2. )
             == Approx( resonances1.front().widths()[2].value ) );
      CHECK( std::sqrt( 1.971943e-2 / 2. )
             == Approx( resonances1.back().widths()[2].value ) );
      CHECK( std::sqrt( 6.012769e-2 / 2. ) == Approx( resonances1.front().eliminatedWidth().value ) );
      CHECK( std::sqrt( 3.926596e-2 / 2. ) == Approx( resonances1.back().eliminatedWidth().value ) );
      */
    } // THEN

    THEN( "cross sections can be reconstructed" ) {

      // values taken from NJOY2016 PeNDF tape for ENDF/B-VII.0 Cf252

      ReactionID elas( "n,Cf252->n,Cf252" );
      ReactionID fiss( "n,Cf252->fission" );
      ReactionID capt( "n,Cf252->capture" );
      std::map< ReactionID, CrossSection > xs;

      xs = resonances( 1e-5 * electronVolt );
      CHECK( 3 == xs.size() );
//      CHECK( 8.152130 == Approx( xs[ elas ].value ) );
//      CHECK( 3.456462e+4 == Approx( xs[ fiss ].value ) );
//      CHECK( 1.284211e+4 == Approx( xs[ capt ].value ) );
    } // THEN
  } // GIVEN
} // SCENARIO

std::string Cf252() {

  // Cf252 ENDF/B-VII.0 LRF=1 resonance evaluation

  return
    " 9.825200+4 2.499160+2          0          0          1          09861 2151     \n"
    " 9.825200+4 1.000000+0          0          1          1          09861 2151     \n"
    " 1.000000-5 3.665000+2          1          1          0          09861 2151     \n"
    " 0.000000+0 8.830960-1          0          0          1          09861 2151     \n"
    " 2.499160+2 0.000000+0          0          0        126         219861 2151     \n"
    "-2.399000+0 5.000000-1 6.164000-2 1.862000-3 2.354000-2 3.624000-29861 2151     \n"
    " 1.711000+1 5.000000-1 8.805000-2 4.972000-3 2.354000-2 5.953000-29861 2151     \n"
    " 3.458000+1 5.000000-1 9.014000-2 7.068000-3 2.354000-2 5.953000-29861 2151     \n"
    " 5.205000+1 5.000000-1 9.174000-2 8.671000-3 2.354000-2 5.953000-29861 2151     \n"
    " 6.952000+1 5.000000-1 9.310000-2 1.002000-2 2.354000-2 5.953000-29861 2151     \n"
    " 8.699000+1 5.000000-1 9.428000-2 1.121000-2 2.354000-2 5.953000-29861 2151     \n"
    " 1.045000+2 5.000000-1 9.536000-2 1.228000-2 2.354000-2 5.953000-29861 2151     \n"
    " 1.219000+2 5.000000-1 9.635000-2 1.327000-2 2.354000-2 5.953000-29861 2151     \n"
    " 1.394000+2 5.000000-1 9.726000-2 1.419000-2 2.354000-2 5.953000-29861 2151     \n"
    " 1.569000+2 5.000000-1 9.813000-2 1.505000-2 2.354000-2 5.953000-29861 2151     \n"
    " 1.743000+2 5.000000-1 9.894000-2 1.587000-2 2.354000-2 5.953000-29861 2151     \n"
    " 1.918000+2 5.000000-1 9.972000-2 1.665000-2 2.354000-2 5.953000-29861 2151     \n"
    " 2.093000+2 5.000000-1 1.005000-1 1.739000-2 2.354000-2 5.953000-29861 2151     \n"
    " 2.267000+2 5.000000-1 1.012000-1 1.810000-2 2.354000-2 5.953000-29861 2151     \n"
    " 2.442000+2 5.000000-1 1.019000-1 1.878000-2 2.354000-2 5.953000-29861 2151     \n"
    " 2.617000+2 5.000000-1 1.025000-1 1.944000-2 2.354000-2 5.953000-29861 2151     \n"
    " 2.792000+2 5.000000-1 1.032000-1 2.008000-2 2.354000-2 5.953000-29861 2151     \n"
    " 2.966000+2 5.000000-1 1.038000-1 2.070000-2 2.354000-2 5.953000-29861 2151     \n"
    " 3.141000+2 5.000000-1 1.044000-1 2.130000-2 2.354000-2 5.953000-29861 2151     \n"
    " 3.316000+2 5.000000-1 1.050000-1 2.189000-2 2.354000-2 5.953000-29861 2151     \n"
    " 3.490000+2 5.000000-1 1.055000-1 2.246000-2 2.354000-2 5.953000-29861 2151     \n"
    " 0.000000+0 0.000000+0          0          0          0          09861 2  0     \n";
}
