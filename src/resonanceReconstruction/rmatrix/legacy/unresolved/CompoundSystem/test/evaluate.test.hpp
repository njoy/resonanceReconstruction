SCENARIO( "evaluate" ) {

  //! @todo add a test with competition

  GIVEN( "Na22 unresolved resonance data" ) {

    // test based on Na22 ENDF/B-VIII.0 LRU=2 unresolved resonance evaluation
    // cross section values extracted from NJOY2016.57

    // this test is essentially for energy independent parameters

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

  GIVEN( "Pu239 unresolved resonance data" ) {

    // test based on Pu239 ENDF/B-VIII.0 LRU=2 unresolved resonance evaluation
    // cross section values extracted from NJOY2016.57

    // this test is essentially for full energy dependent parameters without
    // competition

    // scattering radius from equation D.14
    double a = 0.123 * std::pow( 236.9986 * 1.008664, 1. / 3. ) + 0.08;

    // particles
    Particle neutron( ParticleID( "n" ), neutronMass, 0.0 * coulombs, 0.5, +1);
    Particle pu239( ParticleID( "Pu239" ), 236.9986 * neutronMass,
                    94.0 * coulombs, 0.5, +1); // unknown parity

    // particle pairs
    ParticlePair in( neutron, pu239 );

    // channels
    Channel< Neutron > elastic00( in, in, 0. * electronVolt, { 0, 0.5, 0.0, +1 }, { a * rootBarn, 0.946 * rootBarn } );
    Channel< Neutron > elastic01( in, in, 0. * electronVolt, { 0, 0.5, 1.0, +1 }, { a * rootBarn, 0.946 * rootBarn } );
    Channel< Neutron > elastic10( in, in, 0. * electronVolt, { 1, 0.5, 0.0, +1 }, { a * rootBarn, 0.946 * rootBarn } );
    Channel< Neutron > elastic11( in, in, 0. * electronVolt, { 1, 0.5, 1.0, +1 }, { a * rootBarn, 0.946 * rootBarn } );
    Channel< Neutron > elastic12( in, in, 0. * electronVolt, { 1, 0.5, 2.0, +1 }, { a * rootBarn, 0.946 * rootBarn } );

    // corresponding unresolved resonance tables
    ResonanceTable table00(
      { Resonance( 2.500000e+3 * electronVolt, 8.917200e+0 * electronVolt, 9.508100e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.842000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.550000e+3 * electronVolt, 8.915500e+0 * electronVolt, 8.673000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 4.020000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.650000e+3 * electronVolt, 8.913900e+0 * electronVolt, 1.113200e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.841000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.750000e+3 * electronVolt, 8.912200e+0 * electronVolt, 8.952900e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.840000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.850000e+3 * electronVolt, 8.910500e+0 * electronVolt, 9.309800e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.840000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.950000e+3 * electronVolt, 8.908800e+0 * electronVolt, 1.450600e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 7.840000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.050000e+3 * electronVolt, 8.907200e+0 * electronVolt, 7.845300e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.839000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 3.150000e+3 * electronVolt, 8.905500e+0 * electronVolt, 1.079100e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.838000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 3.250000e+3 * electronVolt, 8.903800e+0 * electronVolt, 1.065500e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.838000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 3.350000e+3 * electronVolt, 8.902200e+0 * electronVolt, 6.778100e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 1.150000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 3.450000e+3 * electronVolt, 8.900500e+0 * electronVolt, 9.239600e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 9.960000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.550000e+3 * electronVolt, 8.898800e+0 * electronVolt, 7.025700e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.836000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 3.650000e+3 * electronVolt, 8.897100e+0 * electronVolt, 1.215900e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 3.930000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.750000e+3 * electronVolt, 8.895500e+0 * electronVolt, 7.980200e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.835000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 3.850000e+3 * electronVolt, 8.893800e+0 * electronVolt, 1.082000e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.835000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 3.950000e+3 * electronVolt, 8.892100e+0 * electronVolt, 1.093700e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.834000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 4.125000e+3 * electronVolt, 8.889100e+0 * electronVolt, 9.226000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.829000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 4.375000e+3 * electronVolt, 8.884900e+0 * electronVolt, 9.985000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.828000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 4.625000e+3 * electronVolt, 8.880800e+0 * electronVolt, 9.217000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.827000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 4.875000e+3 * electronVolt, 8.876600e+0 * electronVolt, 9.974000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.826000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 5.125000e+3 * electronVolt, 8.872400e+0 * electronVolt, 1.011400e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.820000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 5.375000e+3 * electronVolt, 8.868200e+0 * electronVolt, 1.011000e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.820000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 5.625000e+3 * electronVolt, 8.864100e+0 * electronVolt, 1.010500e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.820000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 5.875000e+3 * electronVolt, 8.859900e+0 * electronVolt, 1.010000e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.820000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 6.125000e+3 * electronVolt, 8.855700e+0 * electronVolt, 1.009500e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.820000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 6.375000e+3 * electronVolt, 8.851500e+0 * electronVolt, 1.009100e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.820000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 6.625000e+3 * electronVolt, 8.847300e+0 * electronVolt, 1.008600e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 9.320000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 6.875000e+3 * electronVolt, 8.843100e+0 * electronVolt, 1.008100e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 1.056500e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 7.125000e+3 * electronVolt, 8.839000e+0 * electronVolt, 1.019800e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 4.814000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 7.375000e+3 * electronVolt, 8.834800e+0 * electronVolt, 1.014200e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 1.130000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 7.625000e+3 * electronVolt, 8.830700e+0 * electronVolt, 1.013100e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 1.309500e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 7.875000e+3 * electronVolt, 8.826500e+0 * electronVolt, 1.013400e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.078000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 8.125000e+3 * electronVolt, 8.822300e+0 * electronVolt, 1.005700e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.140000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 8.375000e+3 * electronVolt, 8.818200e+0 * electronVolt, 1.005300e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.807000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 8.625000e+3 * electronVolt, 8.814000e+0 * electronVolt, 1.004800e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.806000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 8.875000e+3 * electronVolt, 8.809900e+0 * electronVolt, 1.004300e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.804000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 9.125000e+3 * electronVolt, 8.805700e+0 * electronVolt, 8.810000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.803000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 9.375000e+3 * electronVolt, 8.801400e+0 * electronVolt, 8.622000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 1.249000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 9.625000e+3 * electronVolt, 8.797400e+0 * electronVolt, 9.643000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 1.250000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 9.875000e+3 * electronVolt, 8.793200e+0 * electronVolt, 9.609000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 1.250000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.025000e+4 * electronVolt, 8.787100e+0 * electronVolt, 9.823000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 1.182000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.075000e+4 * electronVolt, 8.778900e+0 * electronVolt, 9.815000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.794000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.125000e+4 * electronVolt, 8.770500e+0 * electronVolt, 9.805000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.792000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.175000e+4 * electronVolt, 8.762200e+0 * electronVolt, 9.797000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.789000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.225000e+4 * electronVolt, 8.754100e+0 * electronVolt, 9.787000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.786000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.275000e+4 * electronVolt, 8.745800e+0 * electronVolt, 9.778000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.784000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.325000e+4 * electronVolt, 8.737600e+0 * electronVolt, 9.768000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.781000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.375000e+4 * electronVolt, 8.729400e+0 * electronVolt, 9.759000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.779000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.425000e+4 * electronVolt, 8.721200e+0 * electronVolt, 9.750000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 6.620000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.475000e+4 * electronVolt, 8.712900e+0 * electronVolt, 9.740000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.773000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.525000e+4 * electronVolt, 8.704700e+0 * electronVolt, 9.732000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.771000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.575000e+4 * electronVolt, 8.696500e+0 * electronVolt, 9.723000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.768000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.625000e+4 * electronVolt, 8.688300e+0 * electronVolt, 9.714000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.766000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.675000e+4 * electronVolt, 8.680100e+0 * electronVolt, 9.694000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.763000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.725000e+4 * electronVolt, 8.672000e+0 * electronVolt, 9.697000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 9.465000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.775000e+4 * electronVolt, 8.663800e+0 * electronVolt, 9.687000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.758000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.825000e+4 * electronVolt, 8.655700e+0 * electronVolt, 9.678000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.755000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.875000e+4 * electronVolt, 8.647500e+0 * electronVolt, 9.669000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.753000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.925000e+4 * electronVolt, 8.639300e+0 * electronVolt, 9.660000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.750000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.975000e+4 * electronVolt, 8.631200e+0 * electronVolt, 9.651000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.747000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.050000e+4 * electronVolt, 8.619000e+0 * electronVolt, 9.636000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.743000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.150000e+4 * electronVolt, 8.602800e+0 * electronVolt, 9.618000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.738000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.250000e+4 * electronVolt, 8.586700e+0 * electronVolt, 9.600000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.733000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.350000e+4 * electronVolt, 8.570500e+0 * electronVolt, 9.582000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.728000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.450000e+4 * electronVolt, 8.554400e+0 * electronVolt, 9.564000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.723000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.550000e+4 * electronVolt, 8.538200e+0 * electronVolt, 9.546000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.718000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.650000e+4 * electronVolt, 8.522200e+0 * electronVolt, 9.528000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.713000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.750000e+4 * electronVolt, 8.506200e+0 * electronVolt, 9.510000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.708000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.850000e+4 * electronVolt, 8.490100e+0 * electronVolt, 9.492000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.703000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.950000e+4 * electronVolt, 8.474100e+0 * electronVolt, 9.474000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.697000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 3.000000e+4 * electronVolt, 8.465900e+0 * electronVolt, 9.465000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.693000e+0 * electronVolt, 0. * electronVolt )
      },
      { 1, 0, 2, 0 } );
    ResonanceTable table01(
      { Resonance( 2.500000e+3 * electronVolt, 3.044300e+0 * electronVolt, 3.246100e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 7.100000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 2.550000e+3 * electronVolt, 3.043700e+0 * electronVolt, 2.960900e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 3.400000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 2.650000e+3 * electronVolt, 3.043100e+0 * electronVolt, 3.800500e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 1.120000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 2.750000e+3 * electronVolt, 3.042500e+0 * electronVolt, 3.056100e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 1.028000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.850000e+3 * electronVolt, 3.042000e+0 * electronVolt, 3.178300e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 9.160000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 2.950000e+3 * electronVolt, 3.041400e+0 * electronVolt, 4.954000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 1.580000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 3.050000e+3 * electronVolt, 3.040800e+0 * electronVolt, 2.666500e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 5.200000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 3.150000e+3 * electronVolt, 3.040200e+0 * electronVolt, 3.684000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 1.420000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.250000e+3 * electronVolt, 3.039700e+0 * electronVolt, 3.637500e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 6.230000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 3.350000e+3 * electronVolt, 3.039100e+0 * electronVolt, 2.314000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 1.203000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 3.450000e+3 * electronVolt, 3.038500e+0 * electronVolt, 3.154300e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 1.373000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 3.550000e+3 * electronVolt, 3.038000e+0 * electronVolt, 2.398600e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 2.810000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 3.650000e+3 * electronVolt, 3.037400e+0 * electronVolt, 4.151200e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 1.660000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 3.750000e+3 * electronVolt, 3.036800e+0 * electronVolt, 2.724300e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 8.070000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 3.850000e+3 * electronVolt, 3.036200e+0 * electronVolt, 3.693900e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 5.380000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 3.950000e+3 * electronVolt, 3.035700e+0 * electronVolt, 3.733900e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 2.340000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 4.125000e+3 * electronVolt, 3.034600e+0 * electronVolt, 3.148000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 1.098000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 4.375000e+3 * electronVolt, 3.033200e+0 * electronVolt, 3.400000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 2.170000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 4.625000e+3 * electronVolt, 3.031800e+0 * electronVolt, 3.145000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 5.070000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 4.875000e+3 * electronVolt, 3.030400e+0 * electronVolt, 3.403000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 7.500000e-3 * electronVolt, 0. * electronVolt ),
        Resonance( 5.125000e+3 * electronVolt, 3.029000e+0 * electronVolt, 3.453000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 2.360000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 5.375000e+3 * electronVolt, 3.027600e+0 * electronVolt, 3.451000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 1.500000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 5.625000e+3 * electronVolt, 3.026100e+0 * electronVolt, 3.450000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 2.300000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 5.875000e+3 * electronVolt, 3.024700e+0 * electronVolt, 3.448000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 2.220000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 6.125000e+3 * electronVolt, 3.023200e+0 * electronVolt, 3.446000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 3.910000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 6.375000e+3 * electronVolt, 3.021800e+0 * electronVolt, 3.445000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 1.270000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 6.625000e+3 * electronVolt, 3.020400e+0 * electronVolt, 3.443400e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 2.210000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 6.875000e+3 * electronVolt, 3.018900e+0 * electronVolt, 3.441700e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 1.555000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 7.125000e+3 * electronVolt, 3.017500e+0 * electronVolt, 3.483000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 2.920000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 7.375000e+3 * electronVolt, 3.016100e+0 * electronVolt, 3.464000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 6.930000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 7.625000e+3 * electronVolt, 3.014700e+0 * electronVolt, 3.463000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 7.648000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 7.875000e+3 * electronVolt, 3.013300e+0 * electronVolt, 3.463000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 1.544000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 8.125000e+3 * electronVolt, 3.011900e+0 * electronVolt, 3.433000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 4.330000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 8.375000e+3 * electronVolt, 3.010500e+0 * electronVolt, 3.432000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 6.650000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 8.625000e+3 * electronVolt, 3.009000e+0 * electronVolt, 3.430000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 3.180000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 8.875000e+3 * electronVolt, 3.007600e+0 * electronVolt, 3.429000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 5.860000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 9.125000e+3 * electronVolt, 3.006200e+0 * electronVolt, 2.986000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 3.264000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 9.375000e+3 * electronVolt, 3.004800e+0 * electronVolt, 2.923000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 3.136000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 9.625000e+3 * electronVolt, 3.003300e+0 * electronVolt, 3.292000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 7.470000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 9.875000e+3 * electronVolt, 3.001900e+0 * electronVolt, 3.260600e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 3.570000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 1.025000e+4 * electronVolt, 2.999800e+0 * electronVolt, 3.354000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 3.720000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 1.075000e+4 * electronVolt, 2.997000e+0 * electronVolt, 3.351000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 4.610000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 1.125000e+4 * electronVolt, 2.994100e+0 * electronVolt, 3.347000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 3.510000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 1.175000e+4 * electronVolt, 2.991300e+0 * electronVolt, 3.344000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 3.501100e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 1.225000e+4 * electronVolt, 2.988500e+0 * electronVolt, 3.341000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 4.810000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 1.275000e+4 * electronVolt, 2.985700e+0 * electronVolt, 3.337000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 4.710000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 1.325000e+4 * electronVolt, 2.982900e+0 * electronVolt, 3.335000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 4.960000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 1.375000e+4 * electronVolt, 2.980100e+0 * electronVolt, 3.332000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 3.510000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 1.425000e+4 * electronVolt, 2.977300e+0 * electronVolt, 3.329000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 3.660000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 1.475000e+4 * electronVolt, 2.974400e+0 * electronVolt, 3.325000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 5.090000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 1.525000e+4 * electronVolt, 2.971600e+0 * electronVolt, 3.322000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 6.860000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 1.575000e+4 * electronVolt, 2.968800e+0 * electronVolt, 3.319000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 4.220000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 1.625000e+4 * electronVolt, 2.966000e+0 * electronVolt, 3.315000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 6.110000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 1.675000e+4 * electronVolt, 2.963200e+0 * electronVolt, 3.313000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 3.740000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 1.725000e+4 * electronVolt, 2.960400e+0 * electronVolt, 3.310000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 4.130000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 1.775000e+4 * electronVolt, 2.957600e+0 * electronVolt, 3.307000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 8.470000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 1.825000e+4 * electronVolt, 2.954900e+0 * electronVolt, 3.303000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 5.690000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 1.875000e+4 * electronVolt, 2.952100e+0 * electronVolt, 3.301000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 4.590000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 1.925000e+4 * electronVolt, 2.949300e+0 * electronVolt, 3.298000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 6.770000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 1.975000e+4 * electronVolt, 2.946500e+0 * electronVolt, 3.295000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 7.190000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 2.050000e+4 * electronVolt, 2.942400e+0 * electronVolt, 3.290000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 6.100000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 2.150000e+4 * electronVolt, 2.936900e+0 * electronVolt, 3.284000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 6.040000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 2.250000e+4 * electronVolt, 2.931300e+0 * electronVolt, 3.278000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 3.320000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 2.350000e+4 * electronVolt, 2.925800e+0 * electronVolt, 3.271000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 6.640000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 2.450000e+4 * electronVolt, 2.920200e+0 * electronVolt, 3.265000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 7.190000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 2.550000e+4 * electronVolt, 2.914700e+0 * electronVolt, 3.259000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 5.500000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 2.650000e+4 * electronVolt, 2.909300e+0 * electronVolt, 3.253000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 7.870000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 2.750000e+4 * electronVolt, 2.903800e+0 * electronVolt, 3.246000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 6.090000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 2.850000e+4 * electronVolt, 2.898400e+0 * electronVolt, 3.240000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 7.070000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 2.950000e+4 * electronVolt, 2.892900e+0 * electronVolt, 3.234000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 8.710000e-2 * electronVolt, 0. * electronVolt ),
        Resonance( 3.000000e+4 * electronVolt, 2.890100e+0 * electronVolt, 3.231000e-4 * rootElectronVolt, 4.030000e-2 * electronVolt, 7.647000e-2 * electronVolt, 0. * electronVolt )
      },
      { 1, 0, 1, 0 } );
    ResonanceTable table10(
      {
        Resonance( 2.500000e+3 * electronVolt, 8.917200e+0 * electronVolt, 1.573800e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.550000e+3 * electronVolt, 8.915500e+0 * electronVolt, 1.449800e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.650000e+3 * electronVolt, 8.913900e+0 * electronVolt, 1.824400e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.750000e+3 * electronVolt, 8.912200e+0 * electronVolt, 1.497800e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.850000e+3 * electronVolt, 8.910500e+0 * electronVolt, 1.555000e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.950000e+3 * electronVolt, 8.908800e+0 * electronVolt, 2.338500e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 3.050000e+3 * electronVolt, 8.907200e+0 * electronVolt, 1.337900e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 3.150000e+3 * electronVolt, 8.905500e+0 * electronVolt, 1.785600e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 3.250000e+3 * electronVolt, 8.903800e+0 * electronVolt, 1.767600e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 3.350000e+3 * electronVolt, 8.902200e+0 * electronVolt, 1.182400e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 3.450000e+3 * electronVolt, 8.900500e+0 * electronVolt, 1.558300e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 3.550000e+3 * electronVolt, 8.898800e+0 * electronVolt, 1.224400e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 3.650000e+3 * electronVolt, 8.897100e+0 * electronVolt, 2.003600e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 3.750000e+3 * electronVolt, 8.895500e+0 * electronVolt, 1.373800e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 3.850000e+3 * electronVolt, 8.893800e+0 * electronVolt, 1.806300e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 3.950000e+3 * electronVolt, 8.892100e+0 * electronVolt, 1.826000e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 4.125000e+3 * electronVolt, 8.889100e+0 * electronVolt, 1.349000e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 4.375000e+3 * electronVolt, 8.884900e+0 * electronVolt, 1.458600e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 4.625000e+3 * electronVolt, 8.880800e+0 * electronVolt, 1.348000e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 4.875000e+3 * electronVolt, 8.876600e+0 * electronVolt, 1.459600e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 5.125000e+3 * electronVolt, 8.872400e+0 * electronVolt, 1.508300e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 5.375000e+3 * electronVolt, 8.868200e+0 * electronVolt, 1.507600e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 5.625000e+3 * electronVolt, 8.864100e+0 * electronVolt, 1.506900e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 5.875000e+3 * electronVolt, 8.859900e+0 * electronVolt, 1.506200e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 6.125000e+3 * electronVolt, 8.855700e+0 * electronVolt, 1.505500e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 6.375000e+3 * electronVolt, 8.851500e+0 * electronVolt, 1.504800e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 6.625000e+3 * electronVolt, 8.847300e+0 * electronVolt, 1.437500e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 6.875000e+3 * electronVolt, 8.843100e+0 * electronVolt, 1.437500e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 7.125000e+3 * electronVolt, 8.839000e+0 * electronVolt, 1.484100e-3 * rootElectronVolt, 1.170000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 7.375000e+3 * electronVolt, 8.834800e+0 * electronVolt, 1.483400e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 7.625000e+3 * electronVolt, 8.830700e+0 * electronVolt, 1.482800e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 7.875000e+3 * electronVolt, 8.826500e+0 * electronVolt, 1.481100e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 8.125000e+3 * electronVolt, 8.822300e+0 * electronVolt, 1.499800e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 8.375000e+3 * electronVolt, 8.818200e+0 * electronVolt, 1.499100e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 8.625000e+3 * electronVolt, 8.814000e+0 * electronVolt, 1.498400e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 8.875000e+3 * electronVolt, 8.809900e+0 * electronVolt, 1.497700e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 9.125000e+3 * electronVolt, 8.805700e+0 * electronVolt, 1.479500e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 9.375000e+3 * electronVolt, 8.801400e+0 * electronVolt, 1.466200e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 9.625000e+3 * electronVolt, 8.797400e+0 * electronVolt, 1.438700e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 9.875000e+3 * electronVolt, 8.793200e+0 * electronVolt, 1.424500e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.025000e+4 * electronVolt, 8.787100e+0 * electronVolt, 1.463500e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.075000e+4 * electronVolt, 8.778900e+0 * electronVolt, 1.462300e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.125000e+4 * electronVolt, 8.770500e+0 * electronVolt, 1.461100e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.175000e+4 * electronVolt, 8.762200e+0 * electronVolt, 1.459800e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.225000e+4 * electronVolt, 8.754100e+0 * electronVolt, 1.458600e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.275000e+4 * electronVolt, 8.745800e+0 * electronVolt, 1.457400e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.325000e+4 * electronVolt, 8.737600e+0 * electronVolt, 1.456200e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.375000e+4 * electronVolt, 8.729400e+0 * electronVolt, 1.455000e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.425000e+4 * electronVolt, 8.721200e+0 * electronVolt, 1.453700e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.475000e+4 * electronVolt, 8.712900e+0 * electronVolt, 1.452400e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.525000e+4 * electronVolt, 8.704700e+0 * electronVolt, 1.451200e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.575000e+4 * electronVolt, 8.696500e+0 * electronVolt, 1.449900e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.625000e+4 * electronVolt, 8.688300e+0 * electronVolt, 1.448500e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.675000e+4 * electronVolt, 8.680100e+0 * electronVolt, 1.447400e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.725000e+4 * electronVolt, 8.672000e+0 * electronVolt, 1.445700e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.775000e+4 * electronVolt, 8.663800e+0 * electronVolt, 1.444400e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.825000e+4 * electronVolt, 8.655700e+0 * electronVolt, 1.443000e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.875000e+4 * electronVolt, 8.647500e+0 * electronVolt, 1.441700e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.925000e+4 * electronVolt, 8.639300e+0 * electronVolt, 1.440200e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.975000e+4 * electronVolt, 8.631200e+0 * electronVolt, 1.438900e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.050000e+4 * electronVolt, 8.619000e+0 * electronVolt, 1.436800e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.150000e+4 * electronVolt, 8.602800e+0 * electronVolt, 1.434100e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.250000e+4 * electronVolt, 8.586700e+0 * electronVolt, 1.431400e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.350000e+4 * electronVolt, 8.570500e+0 * electronVolt, 1.428700e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.450000e+4 * electronVolt, 8.554400e+0 * electronVolt, 1.426000e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.550000e+4 * electronVolt, 8.538200e+0 * electronVolt, 1.423300e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.650000e+4 * electronVolt, 8.522200e+0 * electronVolt, 1.420600e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.750000e+4 * electronVolt, 8.506200e+0 * electronVolt, 1.418000e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.850000e+4 * electronVolt, 8.490100e+0 * electronVolt, 1.415300e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.950000e+4 * electronVolt, 8.474100e+0 * electronVolt, 1.412600e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 3.000000e+4 * electronVolt, 8.465900e+0 * electronVolt, 1.411300e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt )
      },
      { 1, 0, 0, 0 } );
    ResonanceTable table11(
      {
        Resonance( 2.500000e+3 * electronVolt, 3.044300e+0 * electronVolt, 2.686450e-4 * rootElectronVolt, 3.035000e-2 * electronVolt, 9.720000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.550000e+3 * electronVolt, 3.043700e+0 * electronVolt, 2.474850e-4 * rootElectronVolt, 3.035000e-2 * electronVolt, 9.720000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.650000e+3 * electronVolt, 3.043100e+0 * electronVolt, 3.114150e-4 * rootElectronVolt, 3.035000e-2 * electronVolt, 9.720000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.750000e+3 * electronVolt, 3.042500e+0 * electronVolt, 2.556600e-4 * rootElectronVolt, 3.035000e-2 * electronVolt, 9.710000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.850000e+3 * electronVolt, 3.042000e+0 * electronVolt, 2.654300e-4 * rootElectronVolt, 3.035000e-2 * electronVolt, 9.710000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.950000e+3 * electronVolt, 3.041400e+0 * electronVolt, 3.991700e-4 * rootElectronVolt, 3.035000e-2 * electronVolt, 9.710000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.050000e+3 * electronVolt, 3.040800e+0 * electronVolt, 2.283650e-4 * rootElectronVolt, 3.035000e-2 * electronVolt, 9.710000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.150000e+3 * electronVolt, 3.040200e+0 * electronVolt, 3.047950e-4 * rootElectronVolt, 3.035000e-2 * electronVolt, 9.710000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.250000e+3 * electronVolt, 3.039700e+0 * electronVolt, 3.017200e-4 * rootElectronVolt, 3.035000e-2 * electronVolt, 9.700000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.350000e+3 * electronVolt, 3.039100e+0 * electronVolt, 2.018250e-4 * rootElectronVolt, 3.035000e-2 * electronVolt, 9.700000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.450000e+3 * electronVolt, 3.038500e+0 * electronVolt, 2.659900e-4 * rootElectronVolt, 3.035000e-2 * electronVolt, 9.700000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.550000e+3 * electronVolt, 3.038000e+0 * electronVolt, 2.090000e-4 * rootElectronVolt, 3.035000e-2 * electronVolt, 9.700000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.650000e+3 * electronVolt, 3.037400e+0 * electronVolt, 3.420100e-4 * rootElectronVolt, 3.035000e-2 * electronVolt, 9.700000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.750000e+3 * electronVolt, 3.036800e+0 * electronVolt, 2.345000e-4 * rootElectronVolt, 3.035000e-2 * electronVolt, 9.700000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.850000e+3 * electronVolt, 3.036200e+0 * electronVolt, 3.083250e-4 * rootElectronVolt, 3.035000e-2 * electronVolt, 9.700000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.950000e+3 * electronVolt, 3.035700e+0 * electronVolt, 3.116900e-4 * rootElectronVolt, 3.035000e-2 * electronVolt, 9.700000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 4.125000e+3 * electronVolt, 3.034600e+0 * electronVolt, 2.302000e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.660000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 4.375000e+3 * electronVolt, 3.033200e+0 * electronVolt, 2.488800e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.650000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 4.625000e+3 * electronVolt, 3.031800e+0 * electronVolt, 2.300000e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.650000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 4.875000e+3 * electronVolt, 3.030400e+0 * electronVolt, 2.490000e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.640000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 5.125000e+3 * electronVolt, 3.029000e+0 * electronVolt, 2.574500e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.640000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 5.375000e+3 * electronVolt, 3.027600e+0 * electronVolt, 2.573500e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.640000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 5.625000e+3 * electronVolt, 3.026100e+0 * electronVolt, 2.572000e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.630000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 5.875000e+3 * electronVolt, 3.024700e+0 * electronVolt, 2.571000e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.630000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 6.125000e+3 * electronVolt, 3.023200e+0 * electronVolt, 2.570000e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.620000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 6.375000e+3 * electronVolt, 3.021800e+0 * electronVolt, 2.568500e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.620000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 6.625000e+3 * electronVolt, 3.020400e+0 * electronVolt, 2.452800e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.610000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 6.875000e+3 * electronVolt, 3.018900e+0 * electronVolt, 2.452800e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.610000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 7.125000e+3 * electronVolt, 3.017500e+0 * electronVolt, 2.532500e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.600000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 7.375000e+3 * electronVolt, 3.016100e+0 * electronVolt, 2.531000e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.600000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 7.625000e+3 * electronVolt, 3.014700e+0 * electronVolt, 2.529500e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.600000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 7.875000e+3 * electronVolt, 3.013300e+0 * electronVolt, 2.528500e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.600000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 8.125000e+3 * electronVolt, 3.011900e+0 * electronVolt, 2.560000e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.590000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 8.375000e+3 * electronVolt, 3.010500e+0 * electronVolt, 2.559000e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.580000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 8.625000e+3 * electronVolt, 3.009000e+0 * electronVolt, 2.557500e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.580000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 8.875000e+3 * electronVolt, 3.007600e+0 * electronVolt, 2.556500e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.570000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 9.125000e+3 * electronVolt, 3.006200e+0 * electronVolt, 2.472250e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.570000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 9.375000e+3 * electronVolt, 3.004800e+0 * electronVolt, 2.502650e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.560000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 9.625000e+3 * electronVolt, 3.003300e+0 * electronVolt, 2.455650e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.560000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 9.875000e+3 * electronVolt, 3.001900e+0 * electronVolt, 2.431600e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.550000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.025000e+4 * electronVolt, 2.999800e+0 * electronVolt, 2.500600e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.550000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.075000e+4 * electronVolt, 2.997000e+0 * electronVolt, 2.498200e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.540000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.125000e+4 * electronVolt, 2.994100e+0 * electronVolt, 2.495900e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.530000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.175000e+4 * electronVolt, 2.991300e+0 * electronVolt, 2.493550e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.520000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.225000e+4 * electronVolt, 2.988500e+0 * electronVolt, 2.491150e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.510000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.275000e+4 * electronVolt, 2.985700e+0 * electronVolt, 2.488850e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.500000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.325000e+4 * electronVolt, 2.982900e+0 * electronVolt, 2.486500e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.490000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.375000e+4 * electronVolt, 2.980100e+0 * electronVolt, 2.484200e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.480000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.425000e+4 * electronVolt, 2.977300e+0 * electronVolt, 2.481750e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.480000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.475000e+4 * electronVolt, 2.974400e+0 * electronVolt, 2.479450e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.470000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.525000e+4 * electronVolt, 2.971600e+0 * electronVolt, 2.477100e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.460000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.575000e+4 * electronVolt, 2.968800e+0 * electronVolt, 2.474800e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.450000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.625000e+4 * electronVolt, 2.966000e+0 * electronVolt, 2.472400e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.440000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.675000e+4 * electronVolt, 2.963200e+0 * electronVolt, 2.470100e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.430000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.725000e+4 * electronVolt, 2.960400e+0 * electronVolt, 2.467800e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.420000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.775000e+4 * electronVolt, 2.957600e+0 * electronVolt, 2.465450e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.410000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.825000e+4 * electronVolt, 2.954900e+0 * electronVolt, 2.463100e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.410000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.875000e+4 * electronVolt, 2.952100e+0 * electronVolt, 2.460800e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.400000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.925000e+4 * electronVolt, 2.949300e+0 * electronVolt, 2.458450e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.390000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.975000e+4 * electronVolt, 2.946500e+0 * electronVolt, 2.456150e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.380000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.050000e+4 * electronVolt, 2.942400e+0 * electronVolt, 2.452500e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.360000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.150000e+4 * electronVolt, 2.936900e+0 * electronVolt, 2.448000e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.350000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.250000e+4 * electronVolt, 2.931300e+0 * electronVolt, 2.443500e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.330000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.350000e+4 * electronVolt, 2.925800e+0 * electronVolt, 2.438500e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.310000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.450000e+4 * electronVolt, 2.920200e+0 * electronVolt, 2.434000e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.290000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.550000e+4 * electronVolt, 2.914700e+0 * electronVolt, 2.429500e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.280000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.650000e+4 * electronVolt, 2.904300e+0 * electronVolt, 2.425000e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.260000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.750000e+4 * electronVolt, 2.903800e+0 * electronVolt, 2.420500e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.240000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.850000e+4 * electronVolt, 2.898400e+0 * electronVolt, 2.416000e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.220000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.950000e+4 * electronVolt, 2.892900e+0 * electronVolt, 2.411500e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.210000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.000000e+4 * electronVolt, 2.890100e+0 * electronVolt, 2.408900e-4 * rootElectronVolt, 3.030000e-2 * electronVolt, 9.180000e-1 * electronVolt, 0. * electronVolt )
      },
      { 2, 0, 2, 0 } );
    ResonanceTable table12(
      {
        Resonance( 2.500000e+3 * electronVolt, 1.915900e+0 * electronVolt, 3.381400e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.130000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.550000e+3 * electronVolt, 1.915600e+0 * electronVolt, 3.115100e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.130000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.650000e+3 * electronVolt, 1.915200e+0 * electronVolt, 3.919800e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.120000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.750000e+3 * electronVolt, 1.914800e+0 * electronVolt, 3.218000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.120000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.850000e+3 * electronVolt, 1.914500e+0 * electronVolt, 3.341000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.120000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.950000e+3 * electronVolt, 1.914100e+0 * electronVolt, 5.024300e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.120000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.050000e+3 * electronVolt, 1.913800e+0 * electronVolt, 2.874500e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.120000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.150000e+3 * electronVolt, 1.913400e+0 * electronVolt, 3.836600e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.120000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.250000e+3 * electronVolt, 1.913000e+0 * electronVolt, 3.797700e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.120000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.350000e+3 * electronVolt, 1.912700e+0 * electronVolt, 2.540400e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.120000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.450000e+3 * electronVolt, 1.912200e+0 * electronVolt, 3.347900e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.120000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.550000e+3 * electronVolt, 1.911900e+0 * electronVolt, 2.630600e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.110000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.650000e+3 * electronVolt, 1.911600e+0 * electronVolt, 4.304900e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.110000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.750000e+3 * electronVolt, 1.911200e+0 * electronVolt, 2.951700e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.110000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.850000e+3 * electronVolt, 1.910900e+0 * electronVolt, 3.881000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.110000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.950000e+3 * electronVolt, 1.910500e+0 * electronVolt, 3.923200e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.110000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 4.125000e+3 * electronVolt, 1.909900e+0 * electronVolt, 2.899000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.080000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 4.375000e+3 * electronVolt, 1.909000e+0 * electronVolt, 3.127200e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.080000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 4.625000e+3 * electronVolt, 1.908100e+0 * electronVolt, 2.893000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.070000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 4.875000e+3 * electronVolt, 1.907200e+0 * electronVolt, 3.134000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.070000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 5.125000e+3 * electronVolt, 1.906300e+0 * electronVolt, 3.241000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.070000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 5.375000e+3 * electronVolt, 1.905400e+0 * electronVolt, 3.239000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.070000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 5.625000e+3 * electronVolt, 1.904500e+0 * electronVolt, 3.238000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.070000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 5.875000e+3 * electronVolt, 1.903600e+0 * electronVolt, 3.236000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.060000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 6.125000e+3 * electronVolt, 1.902700e+0 * electronVolt, 3.235000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.060000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 6.375000e+3 * electronVolt, 1.901800e+0 * electronVolt, 3.233000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.060000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 6.625000e+3 * electronVolt, 1.900900e+0 * electronVolt, 3.087800e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.050000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 6.875000e+3 * electronVolt, 1.900000e+0 * electronVolt, 3.087800e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.050000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 7.125000e+3 * electronVolt, 1.899100e+0 * electronVolt, 3.187000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.040000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 7.375000e+3 * electronVolt, 1.898200e+0 * electronVolt, 3.186000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.040000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 7.625000e+3 * electronVolt, 1.897300e+0 * electronVolt, 3.184000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.040000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 7.875000e+3 * electronVolt, 1.896400e+0 * electronVolt, 3.183000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.040000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 8.125000e+3 * electronVolt, 1.895500e+0 * electronVolt, 3.222000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.030000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 8.375000e+3 * electronVolt, 1.894600e+0 * electronVolt, 3.221000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.030000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 8.625000e+3 * electronVolt, 1.893700e+0 * electronVolt, 3.219000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.020000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 8.875000e+3 * electronVolt, 1.892800e+0 * electronVolt, 3.218000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.020000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 9.125000e+3 * electronVolt, 1.891900e+0 * electronVolt, 3.181200e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.020000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 9.375000e+3 * electronVolt, 1.891000e+0 * electronVolt, 3.150100e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.020000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 9.625000e+3 * electronVolt, 1.890100e+0 * electronVolt, 3.091200e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.010000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 9.875000e+3 * electronVolt, 1.889200e+0 * electronVolt, 3.060200e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.010000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.025000e+4 * electronVolt, 1.887900e+0 * electronVolt, 3.147400e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.010000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.075000e+4 * electronVolt, 1.886100e+0 * electronVolt, 3.144500e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.000000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.125000e+4 * electronVolt, 1.884300e+0 * electronVolt, 3.141500e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 6.000000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.175000e+4 * electronVolt, 1.882500e+0 * electronVolt, 3.138500e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.990000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.225000e+4 * electronVolt, 1.880800e+0 * electronVolt, 3.135600e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.990000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.275000e+4 * electronVolt, 1.879000e+0 * electronVolt, 3.132600e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.980000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.325000e+4 * electronVolt, 1.877200e+0 * electronVolt, 3.129600e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.970000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.375000e+4 * electronVolt, 1.875400e+0 * electronVolt, 3.126600e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.970000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.425000e+4 * electronVolt, 1.873700e+0 * electronVolt, 3.123700e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.960000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.475000e+4 * electronVolt, 1.871900e+0 * electronVolt, 3.120700e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.960000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.525000e+4 * electronVolt, 1.870100e+0 * electronVolt, 3.117700e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.950000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.575000e+4 * electronVolt, 1.868400e+0 * electronVolt, 3.114800e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.950000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.625000e+4 * electronVolt, 1.866600e+0 * electronVolt, 3.111900e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.940000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.675000e+4 * electronVolt, 1.864800e+0 * electronVolt, 3.108900e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.940000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.725000e+4 * electronVolt, 1.863100e+0 * electronVolt, 3.106100e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.930000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.775000e+4 * electronVolt, 1.861300e+0 * electronVolt, 3.103100e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.920000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.825000e+4 * electronVolt, 1.859600e+0 * electronVolt, 3.100200e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.920000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.875000e+4 * electronVolt, 1.857800e+0 * electronVolt, 3.097300e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.910000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.925000e+4 * electronVolt, 1.856100e+0 * electronVolt, 3.094400e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.910000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.975000e+4 * electronVolt, 1.854300e+0 * electronVolt, 3.091400e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.900000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.050000e+4 * electronVolt, 1.851700e+0 * electronVolt, 3.087000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.890000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.150000e+4 * electronVolt, 1.848200e+0 * electronVolt, 3.081000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.880000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.250000e+4 * electronVolt, 1.844700e+0 * electronVolt, 3.075000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.870000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.350000e+4 * electronVolt, 1.841300e+0 * electronVolt, 3.070000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.860000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.450000e+4 * electronVolt, 1.837800e+0 * electronVolt, 3.064000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.850000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.550000e+4 * electronVolt, 1.834300e+0 * electronVolt, 3.058000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.840000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.650000e+4 * electronVolt, 1.830800e+0 * electronVolt, 3.052000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.830000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.750000e+4 * electronVolt, 1.827400e+0 * electronVolt, 3.046000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.820000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.850000e+4 * electronVolt, 1.823900e+0 * electronVolt, 3.041000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.800000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.950000e+4 * electronVolt, 1.820500e+0 * electronVolt, 3.035000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.790000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.000000e+4 * electronVolt, 1.818700e+0 * electronVolt, 3.032000e-4 * rootElectronVolt, 3.335000e-2 * electronVolt, 5.770000e-1 * electronVolt, 0. * electronVolt )
      },
      { 1, 0, 2, 0 } );

    // spin groups
    SpinGroup group00( std::move( elastic00 ), std::move( table00 ) );
    SpinGroup group01( std::move( elastic01 ), std::move( table01 ) );
    SpinGroup group10( std::move( elastic10 ), std::move( table10 ) );
    SpinGroup group11( std::move( elastic11 ), std::move( table11 ) );
    SpinGroup group12( std::move( elastic12 ), std::move( table12 ) );

    // the compound system
    CompoundSystem system( { group00, group01, group10, group11, group12 } );

    ReactionID elas( "n,Pu239->n,Pu239" );
    ReactionID capt( "n,Pu239->capture" );
    ReactionID fiss( "n,Pu239->fission" );

    THEN( "cross sections can be calculated" ) {

      std::map< ReactionID, CrossSection > xs;
      system.evaluate( 2500. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 13.431867799494993 == Approx( xs[ elas ].value ) );
      CHECK( 2.4246294118305469 == Approx( xs[ capt ].value ) );
      CHECK( 4.2347077829718600 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 2550. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 13.757049331751272 == Approx( xs[ elas ].value ) );
      CHECK( 2.7535588417225401 == Approx( xs[ capt ].value ) );
      CHECK( 2.7250738402813770 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 2650. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 14.785619814411167 == Approx( xs[ elas ].value ) );
      CHECK( 3.4250992815937775 == Approx( xs[ capt ].value ) );
      CHECK( 3.1034768421167902 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 2750. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 13.015052442868628 == Approx( xs[ elas ].value ) );
      CHECK( 2.0094674966235324 == Approx( xs[ capt ].value ) );
      CHECK( 4.1691468530984688 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 2850. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 13.166067715187607 == Approx( xs[ elas ].value ) );
      CHECK( 2.0766921192491146 == Approx( xs[ capt ].value ) );
      CHECK( 4.1257369135080824 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 2950. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 16.630485325284845 == Approx( xs[ elas ].value ) );
      CHECK( 3.7100196889326873 == Approx( xs[ capt ].value ) );
      CHECK( 3.3622150057314930 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 3050. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.828566425344780 == Approx( xs[ elas ].value ) );
      CHECK( 1.9979093154182013 == Approx( xs[ capt ].value ) );
      CHECK( 3.0165662797466308 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 3150. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 13.383826625287547 == Approx( xs[ elas ].value ) );
      CHECK( 1.9338653188386783 == Approx( xs[ capt ].value ) );
      CHECK( 4.8963477437143519 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 3250. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 13.735444943400172 == Approx( xs[ elas ].value ) );
      CHECK( 2.2765782247111428 == Approx( xs[ capt ].value ) );
      CHECK( 3.9544111307709917 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 3350. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.828163098228208 == Approx( xs[ elas ].value ) );
      CHECK( 2.1659848701277475 == Approx( xs[ capt ].value ) );
      CHECK( 1.7100938479382908 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 3450. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 13.817894408227316 == Approx( xs[ elas ].value ) );
      CHECK( 2.5718330146392030 == Approx( xs[ capt ].value ) );
      CHECK( 2.1984392121425937 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 3550. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.646757244603497 == Approx( xs[ elas ].value ) );
      CHECK( 1.8854034302751068 == Approx( xs[ capt ].value ) );
      CHECK( 2.2139612408447187 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 3650. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 15.310909041429509 == Approx( xs[ elas ].value ) );
      CHECK( 2.9480222219153105 == Approx( xs[ capt ].value ) );
      CHECK( 2.3944569670663327 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 3750. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.640575959272745 == Approx( xs[ elas ].value ) );
      CHECK( 1.6237580373380918 == Approx( xs[ capt ].value ) );
      CHECK( 3.0669295520029860 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 3850. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 13.722029715032093 == Approx( xs[ elas ].value ) );
      CHECK( 2.1218992499625551 == Approx( xs[ capt ].value ) );
      CHECK( 3.5565563633157784 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 3950. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 14.059124787355131 == Approx( xs[ elas ].value ) );
      CHECK( 2.3965504607649137 == Approx( xs[ capt ].value ) );
      CHECK( 2.9310589935473237 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 4125. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 13.528961154554985 == Approx( xs[ elas ].value ) );
      CHECK( 2.2704057957425268 == Approx( xs[ capt ].value ) );
      CHECK( 2.1142575574972580 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 4375. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 13.609780185300750 == Approx( xs[ elas ].value ) );
      CHECK( 2.1293843563176686 == Approx( xs[ capt ].value ) );
      CHECK( 2.5096480505990488 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 4625. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 13.052069928804340 == Approx( xs[ elas ].value ) );
      CHECK( 1.7147979400601743 == Approx( xs[ capt ].value ) );
      CHECK( 2.7719514755312127 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 4875. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 13.719677116547597 == Approx( xs[ elas ].value ) );
      CHECK( 2.1863558896560784 == Approx( xs[ capt ].value ) );
      CHECK( 1.9797770726152923 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 5125. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 13.502623926973984 == Approx( xs[ elas ].value ) );
      CHECK( 1.9162545111510034 == Approx( xs[ capt ].value ) );
      CHECK( 2.4065505126237392 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 5375. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 13.562053958643903 == Approx( xs[ elas ].value ) );
      CHECK( 1.9532105636617532 == Approx( xs[ capt ].value ) );
      CHECK( 2.1533363795826519 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 5625. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 13.422695511184223 == Approx( xs[ elas ].value ) );
      CHECK( 1.8068452088555180 == Approx( xs[ capt ].value ) );
      CHECK( 2.2939498019909874 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 5875. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 13.390192639539050 == Approx( xs[ elas ].value ) );
      CHECK( 1.7627114148361509 == Approx( xs[ capt ].value ) );
      CHECK( 2.2332467864787167 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 6125. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 13.198804329211262 == Approx( xs[ elas ].value ) );
      CHECK( 1.5843835161262598 == Approx( xs[ capt ].value ) );
      CHECK( 2.4741616941079800 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 6375. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 13.424706176346824 == Approx( xs[ elas ].value ) );
      CHECK( 1.7678398184607547 == Approx( xs[ capt ].value ) );
      CHECK( 1.9446042404062127 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 6625. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 13.455973570069764 == Approx( xs[ elas ].value ) );
      CHECK( 1.6667307768985120 == Approx( xs[ capt ].value ) );
      CHECK( 1.8811526705591437 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 6875. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 13.467089082469917 == Approx( xs[ elas ].value ) );
      CHECK( 1.6814347745297269 == Approx( xs[ capt ].value ) );
      CHECK( 1.7470805116987502 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 7125. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 13.488548279720980 == Approx( xs[ elas ].value ) );
      CHECK( 1.5784312329223060 == Approx( xs[ capt ].value ) );
      CHECK( 1.8044139743918477 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 7375. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 13.020788627011836 == Approx( xs[ elas ].value ) );
      CHECK( 1.3053812558314954 == Approx( xs[ capt ].value ) );
      CHECK( 2.4205352196525434 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 7625. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.934090129198637 == Approx( xs[ elas ].value ) );
      CHECK( 1.2483261737077900 == Approx( xs[ capt ].value ) );
      CHECK( 2.4706296545479458 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 7875. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 13.246671145851739 == Approx( xs[ elas ].value ) );
      CHECK( 1.5201992860824949 == Approx( xs[ capt ].value ) );
      CHECK( 1.7997565204799648 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 8125. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.958211772802807 == Approx( xs[ elas ].value ) );
      CHECK( 1.3036570594495487 == Approx( xs[ capt ].value ) );
      CHECK( 2.1867221998170163 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 8375. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.770586335936564 == Approx( xs[ elas ].value ) );
      CHECK( 1.1804309950377845 == Approx( xs[ capt ].value ) );
      CHECK( 2.4176316769531470 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 8625. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.947806523947804 == Approx( xs[ elas ].value ) );
      CHECK( 1.3050046097936769 == Approx( xs[ capt ].value ) );
      CHECK( 2.0381641692884105 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate(  8875. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.760432686695797 == Approx( xs[ elas ].value ) );
      CHECK( 1.1640784418020065 == Approx( xs[ capt ].value ) );
      CHECK( 2.2929613239197519 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate(  9125. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.559170184948195 == Approx( xs[ elas ].value ) );
      CHECK( 1.1655220356423068 == Approx( xs[ capt ].value ) );
      CHECK( 1.8347520770291275 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate(  9375. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.592262067702467 == Approx( xs[ elas ].value ) );
      CHECK( 1.1595571193848930 == Approx( xs[ capt ].value ) );
      CHECK( 1.6663681889910804 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate(  9625. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.635026246096960 == Approx( xs[ elas ].value ) );
      CHECK( 1.0544019009565162 == Approx( xs[ capt ].value ) );
      CHECK( 2.1334064846186003 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate(  9875. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.782067401286186 == Approx( xs[ elas ].value ) );
      CHECK( 1.1648992178296613 == Approx( xs[ capt ].value ) );
      CHECK( 1.7784561774099026 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 10250. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.819335366868941 == Approx( xs[ elas ].value ) );
      CHECK( 1.1507794220832457 == Approx( xs[ capt ].value ) );
      CHECK( 1.7951686104929920 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 10750. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.611896798348319 == Approx( xs[ elas ].value ) );
      CHECK( 1.0596920794676248 == Approx( xs[ capt ].value ) );
      CHECK( 1.9818319032543621 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 11250. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.633103421665119 == Approx( xs[ elas ].value ) );
      CHECK( 1.0711604035852116 == Approx( xs[ capt ].value ) );
      CHECK( 1.8434964727682279 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 11750. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.594176815248932 == Approx( xs[ elas ].value ) );
      CHECK( 1.0419874775673770 == Approx( xs[ capt ].value ) );
      CHECK( 1.8133049719717498 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 12250. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.486667608800810 == Approx( xs[ elas ].value ) );
      CHECK( 0.96866065027378934 == Approx( xs[ capt ].value ) );
      CHECK(  1.901355820300749 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 12750. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.455471187853776 == Approx( xs[ elas ].value ) );
      CHECK( 0.94742843284310452 == Approx( xs[ capt ].value ) );
      CHECK(  1.865477512557371 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 13250. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.410143794353386 == Approx( xs[ elas ].value ) );
      CHECK( 0.91755631630504819 == Approx( xs[ capt ].value ) );
      CHECK(  1.858572995399656 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 13750. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.449135602633906 == Approx( xs[ elas ].value ) );
      CHECK( 0.94201061234645433 == Approx( xs[ capt ].value ) );
      CHECK(  1.716542591252183 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 14250. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.591634237352824 == Approx( xs[ elas ].value ) );
      CHECK( 0.94850656944138445 == Approx( xs[ capt ].value ) );
      CHECK(  1.492806266937279 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 14750. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.308022371357769 == Approx( xs[ elas ].value ) );
      CHECK( 0.85450571851375257 == Approx( xs[ capt ].value ) );
      CHECK(  1.798497061940447 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 15250. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.210845244801529 == Approx( xs[ elas ].value ) );
      CHECK( 0.79752517089583075 == Approx( xs[ capt ].value ) );
      CHECK(  1.884810485963839 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 15750. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.286799583842889 == Approx( xs[ elas ].value ) );
      CHECK( 0.84320542472853222 == Approx( xs[ capt ].value ) );
      CHECK(  1.698221819747725 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 16250. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.180623078524269 == Approx( xs[ elas ].value ) );
      CHECK( 0.78202638203132013 == Approx( xs[ capt ].value ) );
      CHECK(  1.802674092760821 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 16750. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.251610953569406 == Approx( xs[ elas ].value ) );
      CHECK( 0.82483223820061735 == Approx( xs[ capt ].value ) );
      CHECK(  1.629080574217000 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 17250. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.329807047967776 == Approx( xs[ elas ].value ) );
      CHECK( 0.81924952671935924 == Approx( xs[ capt ].value ) );
      CHECK(  1.500099191463520 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 17750. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.028767354762378 == Approx( xs[ elas ].value ) );
      CHECK( 0.70155388127063223 == Approx( xs[ capt ].value ) );
      CHECK(  1.863657774693272 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 18250. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.090362225314085 == Approx( xs[ elas ].value ) );
      CHECK( 0.73696423297333669 == Approx( xs[ capt ].value ) );
      CHECK(  1.712835190216758 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 18750. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.107455154053609 == Approx( xs[ elas ].value ) );
      CHECK( 0.74833832603990780 == Approx( xs[ capt ].value ) );
      CHECK(  1.633891650953657 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 19250. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 12.005872458454432 == Approx( xs[ elas ].value ) );
      CHECK( 0.69448668036096795 == Approx( xs[ capt ].value ) );
      CHECK(  1.739938668328656 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 19750. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 11.969798878422685 == Approx( xs[ elas ].value ) );
      CHECK( 0.67732933688504371 == Approx( xs[ capt ].value ) );
      CHECK(  1.745465190754702 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 20500. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 11.968800519801665 == Approx( xs[ elas ].value ) );
      CHECK( 0.68008741025293007 == Approx( xs[ capt ].value ) );
      CHECK(  1.674538629296025 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 21500. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 11.925981702969445 == Approx( xs[ elas ].value ) );
      CHECK( 0.66213720180738678 == Approx( xs[ capt ].value ) );
      CHECK(  1.648644743940811 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 22500. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 11.983238513168740 == Approx( xs[ elas ].value ) );
      CHECK( 0.69807106458315982 == Approx( xs[ capt ].value ) );
      CHECK(  1.473934886435610 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 23500. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 11.823209511276264 == Approx( xs[ elas ].value ) );
      CHECK( 0.61943892610228057 == Approx( xs[ capt ].value ) );
      CHECK(  1.635054817328119 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 24500. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 11.768399134539791 == Approx( xs[ elas ].value ) );
      CHECK( 0.59731025362879375 == Approx( xs[ capt ].value ) );
      CHECK(  1.639098826246428 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 25500. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 11.778121013407796 == Approx( xs[ elas ].value ) );
      CHECK( 0.60740269130167657 == Approx( xs[ capt ].value ) );
      CHECK(  1.550000939633398 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 26500. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 11.676168727767354 == Approx( xs[ elas ].value ) );
      CHECK( 0.56319077036259035 == Approx( xs[ capt ].value ) );
      CHECK(  1.630580402041155 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 27500. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 11.685793766139076 == Approx( xs[ elas ].value ) );
      CHECK( 0.57319551308004868 == Approx( xs[ capt ].value ) );
      CHECK(  1.547042997823886 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 28500. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 11.624645091802503 == Approx( xs[ elas ].value ) );
      CHECK( 0.54965485457342067 == Approx( xs[ capt ].value ) );
      CHECK(  1.571659583301897 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 29500. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 11.553222436781994 == Approx( xs[ elas ].value ) );
      CHECK( 0.52209263866606681 == Approx( xs[ capt ].value ) );
      CHECK(  1.612989558297598 == Approx( xs[ fiss ].value ) );
      xs.clear();

      system.evaluate( 30000. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 11.560334256440548 == Approx( xs[ elas ].value ) );
      CHECK( 0.52805263879204856 == Approx( xs[ capt ].value ) );
      CHECK(  1.572007382214403 == Approx( xs[ fiss ].value ) );
      xs.clear();
    } // THEN
  } // GIVEN
} // SCENARIO
