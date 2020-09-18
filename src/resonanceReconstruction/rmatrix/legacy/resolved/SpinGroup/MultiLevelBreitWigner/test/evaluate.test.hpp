SCENARIO( "evaluate" ) {

  GIVEN( "Rh105 resolved resonance data" ) {

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

    // spin group
    SpinGroup< MultiLevelBreitWigner > group( std::move( elastic ), std::move( table ), 0. * electronVolt );

    ReactionID elas( "n,Rh105->n,Rh105" );
    ReactionID capt( "n,Rh105->capture" );

    THEN( "cross sections can be calculated for l=0" ) {

      std::map< ReactionID, CrossSection > xs;
      group.evaluate( 1e-5 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 5268.5966369331500 == Approx( xs[ elas ].value ) );
      CHECK( 801565.16324338294  == Approx( xs[ capt ].value ) );
      xs.clear();
/*
      group.evaluate( 1e-4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 9070.5097508583694 == Approx( xs[ elas ].value ) );
      CHECK( 253468.26154893031 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e-3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 9067.1505748817472 == Approx( xs[ elas ].value ) );
      CHECK( 80132.232264999941 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e-2 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 9036.3652957521754 == Approx( xs[ elas ].value ) );
      CHECK( 25279.102640176850 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e-1 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8746.8273681445971 == Approx( xs[ elas ].value ) );
      CHECK( 7816.7206941266395 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 6653.1682520698705 == Approx( xs[ elas ].value ) );
      CHECK( 2161.1909561504717 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 4.755 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 93347.351129302318 == Approx( xs[ elas ].value ) );
      CHECK( 45893.526706445802 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 5. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 182825.08293451308 == Approx( xs[ elas ].value ) );
      CHECK( 87782.287793509589 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 5.245 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 91772.055495597771 == Approx( xs[ elas ].value ) );
      CHECK( 42264.081005233311 == Approx( xs[ capt ].value ) );
      xs.clear();*/
    } // THEN
  } // GIVEN
} // SCENARIO
