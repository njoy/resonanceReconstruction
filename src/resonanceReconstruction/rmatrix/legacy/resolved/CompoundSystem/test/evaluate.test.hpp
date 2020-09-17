SCENARIO( "evaluate" ) {

  //! @todo add a test with competition and fission

  GIVEN( "Rh105 resolved resonance data" ) {

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
    } // THEN
  } // GIVEN
} // SCENARIO
