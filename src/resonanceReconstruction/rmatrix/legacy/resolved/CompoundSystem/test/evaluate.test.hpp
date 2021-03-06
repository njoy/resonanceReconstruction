SCENARIO( "evaluate" ) {

  //! @todo add a test with competition and fission

  GIVEN( "Rh105 resolved resonance data using SLBW" ) {

    // test based on Rh105 ENDF/B-VII.1 resolved resonance evaluation
    // cross section values extracted from NJOY2016.57

    // scattering radius from equation D.14
    double a = 0.123 * std::pow( 104.005 * 1.008664, 1. / 3. ) + 0.08;

    // particles
    Particle neutron( ParticleID( "n" ), neutronMass, 0.0 * coulombs, 0.5, +1);
    Particle rh105( ParticleID( "Rh105" ), 104.005 * neutronMass,
                    45.0 * elementary, 0.5, +1);

    // particle pairs
    ParticlePair in( neutron, rh105 );

    // channels
    Channel< Neutron > elastic( in, in, 0. * electronVolt,
                                { 0, 0.5, 1.0, +1 },
                                { a * rootBarn, 0.62 * rootBarn } );

    // resonance table
    ResonanceTable table(
      { Resonance( -5. * electronVolt,
                   1.45 * electronVolt, 0.16 * electronVolt,
                   0. * electronVolt, 0. * electronVolt,
                   elastic.penetrability( -5. * electronVolt ),
                   elastic.penetrability( -5. * electronVolt ),
                   elastic.shiftFactor( -5. * electronVolt ) ),
        Resonance( 5. * electronVolt,
                   0.33 * electronVolt, 0.16 * electronVolt,
                   0. * electronVolt, 0. * electronVolt,
                   elastic.penetrability( 5. * electronVolt ),
                   elastic.penetrability( 5. * electronVolt ),
                   elastic.shiftFactor( 5. * electronVolt ) ) } );

    SpinGroup< SingleLevelBreitWigner > sgroup( std::move( elastic ), std::move( table ), 0. * electronVolt );

    // the compound system
    CompoundSystem< SingleLevelBreitWigner > system( { sgroup } );

    ReactionID elas( "n,Rh105->n,Rh105" );
    ReactionID capt( "n,Rh105->capture" );

    THEN( "cross sections can be calculated" ) {

      std::map< ReactionID, CrossSection > xs;
      system.evaluate( 1e-5 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 9075.762 == Approx( xs[ elas ].value ) );
      CHECK( 801565.16324338294  == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e-4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 9075.340 == Approx( xs[ elas ].value ) );
      CHECK( 253468.26154893031 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e-3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 9071.981 == Approx( xs[ elas ].value ) );
      CHECK( 80132.232264999941 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e-2 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 9041.196 == Approx( xs[ elas ].value ) );
      CHECK( 25279.102640176850 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e-1 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8751.658 == Approx( xs[ elas ].value ) );
      CHECK( 7816.7206941266395 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 6657.999 == Approx( xs[ elas ].value ) );
      CHECK( 2161.1909561504717 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 4.755 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 9.335218e+4 == Approx( xs[ elas ].value ) );
      CHECK( 45893.526706445802 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 5. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 1.828299e+5 == Approx( xs[ elas ].value ) );
      CHECK( 87782.287793509589 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 5.245 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 9.177689e+4 == Approx( xs[ elas ].value ) );
      CHECK( 42264.081005233311 == Approx( xs[ capt ].value ) );
      xs.clear();
    } // THEN
  } // GIVEN

  GIVEN( "Rh105 resolved resonance data using MLBW" ) {

    // test based on Rh105 ENDF/B-VII.1 resolved resonance evaluation
    // cross section values extracted from NJOY2016.57

    // scattering radius from equation D.14
    double a = 0.123 * std::pow( 104.005 * 1.008664, 1. / 3. ) + 0.08;

    // particles
    Particle neutron( ParticleID( "n" ), neutronMass, 0.0 * coulombs, 0.5, +1);
    Particle rh105( ParticleID( "Rh105" ), 104.005 * neutronMass,
                    45.0 * elementary, 0.5, +1);

    // particle pairs
    ParticlePair in( neutron, rh105 );

    // channels
    Channel< Neutron > elastic( in, in, 0. * electronVolt,
                                { 0, 0.5, 1.0, +1 },
                                { a * rootBarn, 0.62 * rootBarn } );

    // resonance table
    ResonanceTable table(
      { Resonance( -5. * electronVolt,
                   1.45 * electronVolt, 0.16 * electronVolt,
                   0. * electronVolt, 0. * electronVolt,
                   elastic.penetrability( -5. * electronVolt ),
                   elastic.penetrability( -5. * electronVolt ),
                   elastic.shiftFactor( -5. * electronVolt ) ),
        Resonance( 5. * electronVolt,
                   0.33 * electronVolt, 0.16 * electronVolt,
                   0. * electronVolt, 0. * electronVolt,
                   elastic.penetrability( 5. * electronVolt ),
                   elastic.penetrability( 5. * electronVolt ),
                   elastic.shiftFactor( 5. * electronVolt ) ) } );

    SpinGroup< MultiLevelBreitWigner > sgroup( std::move( elastic ), std::move( table ), 0. * electronVolt );

    // the compound system
    CompoundSystem< MultiLevelBreitWigner > system( { sgroup } );

    ReactionID elas( "n,Rh105->n,Rh105" );
    ReactionID capt( "n,Rh105->capture" );

    THEN( "cross sections can be calculated" ) {

      std::map< ReactionID, CrossSection > xs;
      system.evaluate( 1e-5 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 5268.5966369331500 == Approx( xs[ elas ].value ) );
      CHECK( 801565.16324338294  == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e-4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 5268.2765170745824 == Approx( xs[ elas ].value ) );
      CHECK( 253468.26154893031 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e-3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 5265.2500341956293 == Approx( xs[ elas ].value ) );
      CHECK( 80132.232264999941 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e-2 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 5235.6252833826566 == Approx( xs[ elas ].value ) );
      CHECK( 25279.102640176850 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e-1 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4949.4375301958371 == Approx( xs[ elas ].value ) );
      CHECK( 7816.7206941266395 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2723.0907713948382 == Approx( xs[ elas ].value ) );
      CHECK( 2161.1909561504717 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 4.755 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 74816.910028973463 == Approx( xs[ elas ].value ) );
      CHECK( 45893.526706445802 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 5. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 185939.51714628452 == Approx( xs[ elas ].value ) );
      CHECK( 87782.287793509589 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 5.245 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 111838.57454107259 == Approx( xs[ elas ].value ) );
      CHECK( 42264.081005233311 == Approx( xs[ capt ].value ) );
      xs.clear();
    } // THEN
  } // GIVEN
} // SCENARIO
