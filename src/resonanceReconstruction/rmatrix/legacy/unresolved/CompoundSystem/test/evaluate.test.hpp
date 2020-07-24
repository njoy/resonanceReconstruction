SCENARIO( "evaluate" ) {

  GIVEN( "Na22 unresolved resonance data" ) {

    // test based on Na22 ENDF/B-VIII.0 LRU=2 unresolved resonance evaluation
    // cross section values extracted from NJOY2016.57

    // scattering radius from equation D.14
    double a = 0.123 * std::pow( 21.80550 * 1.008664, 1. / 3. ) + 0.08;

    // particles
    Particle neutron( ParticleID( "n" ), neutronMass, 0.0 * coulombs, 0.5, +1);
    Particle na22( ParticleID( "Na22" ), 21.80550 * neutronMass,
                   11.0 * coulombs, 3.0, +1); // unknown parity

    // particle pairs
    ParticlePair in( neutron, na22 );

    // channels
    Channel< Neutron > elastic00( in, in, 0. * electronVolt,
                                  { 0, 0.5, 2.5, +1 },
                                  { a * rootBarn, 0.57 * rootBarn } );
    Channel< Neutron > elastic01( in, in, 0. * electronVolt,
                                  { 0, 0.5, 3.5, +1 },
                                  { a * rootBarn, 0.57 * rootBarn } );
    Channel< Neutron > elastic10( in, in, 0. * electronVolt,
                                  { 1, 0.5, 1.5, +1 },
                                  { a * rootBarn, 0.57 * rootBarn } );
    Channel< Neutron > elastic11( in, in, 0. * electronVolt,
                                  { 1, 0.5, 2.5, +1 },
                                  { a * rootBarn, 0.57 * rootBarn } );
    Channel< Neutron > elastic12( in, in, 0. * electronVolt,
                                  { 1, 0.5, 3.5, +1 },
                                  { a * rootBarn, 0.57 * rootBarn } );
    Channel< Neutron > elastic13( in, in, 0. * electronVolt,
                                  { 1, 0.5, 4.5, +1 },
                                  { a * rootBarn, 0.57 * rootBarn } );
    Channel< Neutron > elastic20( in, in, 0. * electronVolt,
                                  { 2, 0.5, 0.5, +1 },
                                  { a * rootBarn, 0.57 * rootBarn } );
    Channel< Neutron > elastic21( in, in, 0. * electronVolt,
                                  { 2, 0.5, 1.5, +1 },
                                  { a * rootBarn, 0.57 * rootBarn } );
    Channel< Neutron > elastic22( in, in, 0. * electronVolt,
                                  { 2, 0.5, 2.5, +1 },
                                  { a * rootBarn, 0.57 * rootBarn } );
    Channel< Neutron > elastic23( in, in, 0. * electronVolt,
                                  { 2, 0.5, 3.5, +1 },
                                  { a * rootBarn, 0.57 * rootBarn } );
    Channel< Neutron > elastic24( in, in, 0. * electronVolt,
                                  { 2, 0.5, 4.5, +1 },
                                  { a * rootBarn, 0.57 * rootBarn } );
    Channel< Neutron > elastic25( in, in, 0. * electronVolt,
                                  { 2, 0.5, 5.5, +1 },
                                  { a * rootBarn, 0.57 * rootBarn } );

    // corresponding unresolved resonance tables
    ResonanceTable table00(
      { Resonance( 1.5e+4 * electronVolt, 55793.70 * electronVolt, 7.476360 * rootElectronVolt, 1.081650 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.0e+5 * electronVolt, 55793.70 * electronVolt, 7.476360 * rootElectronVolt, 1.081650 * electronVolt, 0. * electronVolt, 0. * electronVolt ) },
      { 1, 0, 0, 0 } );
    ResonanceTable table01(
      { Resonance( 1.5e+4 * electronVolt, 68407.20 * electronVolt, 9.166560 * rootElectronVolt, 1.081650 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.0e+5 * electronVolt, 68407.20 * electronVolt, 9.166560 * rootElectronVolt, 1.081650 * electronVolt, 0. * electronVolt, 0. * electronVolt ) },
      { 1, 0, 0, 0 } );
    ResonanceTable table10(
      { Resonance( 1.5e+4 * electronVolt, 58912.60 * electronVolt, 16.43660 * rootElectronVolt, 5.444310 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.0e+5 * electronVolt, 58912.60 * electronVolt, 16.43660 * rootElectronVolt, 5.444310 * electronVolt, 0. * electronVolt, 0. * electronVolt ) },
      { 1, 0, 0, 0 } );
    ResonanceTable table11(
      { Resonance( 1.5e+4 * electronVolt, 55793.70 * electronVolt, 31.13290 * rootElectronVolt, 5.444310 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.0e+5 * electronVolt, 55793.70 * electronVolt, 31.13290 * rootElectronVolt, 5.444310 * electronVolt, 0. * electronVolt, 0. * electronVolt ) },
      { 2, 0, 0, 0 } );
    ResonanceTable table12(
      { Resonance( 1.5e+4 * electronVolt, 68407.20 * electronVolt, 38.17120 * rootElectronVolt, 5.444310 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.0e+5 * electronVolt, 68407.20 * electronVolt, 38.17120 * rootElectronVolt, 5.444310 * electronVolt, 0. * electronVolt, 0. * electronVolt ) },
      { 2, 0, 0, 0 } );
    ResonanceTable table13(
      { Resonance( 1.5e+4 * electronVolt, 102952.0 * electronVolt, 28.72350 * rootElectronVolt, 5.444310 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.0e+5 * electronVolt, 102952.0 * electronVolt, 28.72350 * rootElectronVolt, 5.444310 * electronVolt, 0. * electronVolt, 0. * electronVolt ) },
      { 1, 0, 0, 0 } );
    ResonanceTable table20(
      { Resonance( 1.5e+4 * electronVolt, 95446.00 * electronVolt, 31.49720 * rootElectronVolt, 1.081650 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.0e+5 * electronVolt, 95446.00 * electronVolt, 31.49720 * rootElectronVolt, 1.081650 * electronVolt, 0. * electronVolt, 0. * electronVolt ) },
      { 1, 0, 0, 0 } );
    ResonanceTable table21(
      { Resonance( 1.5e+4 * electronVolt, 58912.60 * electronVolt, 38.88230 * rootElectronVolt, 1.081650 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.0e+5 * electronVolt, 58912.60 * electronVolt, 38.88230 * rootElectronVolt, 1.081650 * electronVolt, 0. * electronVolt, 0. * electronVolt ) },
      { 2, 0, 0, 0 } );
    ResonanceTable table22(
      { Resonance( 1.5e+4 * electronVolt, 55793.70 * electronVolt, 36.82380 * rootElectronVolt, 1.081650 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.0e+5 * electronVolt, 55793.70 * electronVolt, 36.82380 * rootElectronVolt, 1.081650 * electronVolt, 0. * electronVolt, 0. * electronVolt ) },
      { 2, 0, 0, 0 } );
    ResonanceTable table23(
      { Resonance( 1.5e+4 * electronVolt, 68407.20 * electronVolt, 45.14870 * rootElectronVolt, 1.081650 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.0e+5 * electronVolt, 68407.20 * electronVolt, 45.14870 * rootElectronVolt, 1.081650 * electronVolt, 0. * electronVolt, 0. * electronVolt ) },
      { 2, 0, 0, 0 } );
    ResonanceTable table24(
      { Resonance( 1.5e+4 * electronVolt, 102952.0 * electronVolt, 67.94820 * rootElectronVolt, 1.081650 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.0e+5 * electronVolt, 102952.0 * electronVolt, 67.94820 * rootElectronVolt, 1.081650 * electronVolt, 0. * electronVolt, 0. * electronVolt ) },
      { 2, 0, 0, 0 } );
    ResonanceTable table25(
      { Resonance( 1.5e+4 * electronVolt, 185730.0 * electronVolt, 61.29090 * rootElectronVolt, 1.081650 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.0e+5 * electronVolt, 185730.0 * electronVolt, 61.29090 * rootElectronVolt, 1.081650 * electronVolt, 0. * electronVolt, 0. * electronVolt ) },
      { 1, 0, 0, 0 } );

    // spin groups
    SpinGroup group00( std::move( elastic00 ), std::move( table00 ) );
    SpinGroup group01( std::move( elastic01 ), std::move( table01 ) );
    SpinGroup group10( std::move( elastic10 ), std::move( table10 ) );
    SpinGroup group11( std::move( elastic11 ), std::move( table11 ) );
    SpinGroup group12( std::move( elastic12 ), std::move( table12 ) );
    SpinGroup group13( std::move( elastic13 ), std::move( table13 ) );
    SpinGroup group20( std::move( elastic20 ), std::move( table20 ) );
    SpinGroup group21( std::move( elastic21 ), std::move( table21 ) );
    SpinGroup group22( std::move( elastic22 ), std::move( table22 ) );
    SpinGroup group23( std::move( elastic23 ), std::move( table23 ) );
    SpinGroup group24( std::move( elastic24 ), std::move( table24 ) );
    SpinGroup group25( std::move( elastic25 ), std::move( table25 ) );

    // the compound system
    CompoundSystem system( { group00, group01, group10, group11, group12,
                             group13, group20, group21, group22, group23,
                             group24, group25 } );

    ReactionID elas( "n,Na22->n,Na22" );
    ReactionID capt( "n,Na22->capture" );

    THEN( "cross sections can be calculated" ) {

      std::map< ReactionID, CrossSection > xs;
      system.evaluate( 1.5e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 9.3027405025178052 == Approx( xs[ elas ].value ) );
      CHECK( 4.0639184850639387E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 2e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.7051878520691872 == Approx( xs[ elas ].value ) );
      CHECK( 3.2650798234536135E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 3e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.0164891832868452 == Approx( xs[ elas ].value ) );
      CHECK( 2.3797442539692449E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 4e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 7.6215788391323818 == Approx( xs[ elas ].value ) );
      CHECK( 1.8870640640168301E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 5e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 7.3617143134432972 == Approx( xs[ elas ].value ) );
      CHECK( 1.5670838905457892E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 6e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 7.1763623147011355 == Approx( xs[ elas ].value ) );
      CHECK( 1.3406527410576555E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 7e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 7.0369090178438851 == Approx( xs[ elas ].value ) );
      CHECK( 1.1714340063536310E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 8e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 6.9279078516440951 == Approx( xs[ elas ].value ) );
      CHECK( 1.0400317466079910E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 9e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 6.8402239876247961 == Approx( xs[ elas ].value ) );
      CHECK( 9.3501533213065698E-003 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e+5 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 6.7680817547271817 == Approx( xs[ elas ].value ) );
      CHECK( 8.4916589760332239E-003 == Approx( xs[ capt ].value ) );
      xs.clear();
    } // THEN
  } // GIVEN
} // SCENARIO
