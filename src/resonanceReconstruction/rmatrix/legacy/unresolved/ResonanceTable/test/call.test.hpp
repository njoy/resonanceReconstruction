SCENARIO( "call operator" ) {

  GIVEN( "a valid ResonanceTable for energy independent parameters" ) {

    ResonanceTable table( { { 1000. * electronVolt,
                               250. * electronVolt,
                               100. * rootElectronVolt,
                               500. * electronVolt,
                               400. * electronVolt,
                               200. * electronVolt },
                            { 5000. * electronVolt,
                               250. * electronVolt,
                               100. * rootElectronVolt,
                               500. * electronVolt,
                               400. * electronVolt,
                               200. * electronVolt } },
                          { 2, 0, 4, 1 } );

    THEN( "Unresolved resonances can be retrieved at any energy" ) {

      Resonance resonance = table( 1000. * electronVolt );

      CHECK( 1000. == Approx( resonance.energy().value ) );
      CHECK( 250. == Approx( resonance.levelSpacing().value ) );
      CHECK( 100. == Approx( resonance.elastic().value ) );
      CHECK( 500. == Approx( resonance.capture().value ) );
      CHECK( 400. == Approx( resonance.fission().value ) );
      CHECK( 200. == Approx( resonance.competition().value ) );

      resonance = table( 3000. * electronVolt );

      CHECK( 3000. == Approx( resonance.energy().value ) );
      CHECK( 250. == Approx( resonance.levelSpacing().value ) );
      CHECK( 100. == Approx( resonance.elastic().value ) );
      CHECK( 500. == Approx( resonance.capture().value ) );
      CHECK( 400. == Approx( resonance.fission().value ) );
      CHECK( 200. == Approx( resonance.competition().value ) );

      resonance = table( 5000. * electronVolt );

      CHECK( 5000. == Approx( resonance.energy().value ) );
      CHECK( 250. == Approx( resonance.levelSpacing().value ) );
      CHECK( 100. == Approx( resonance.elastic().value ) );
      CHECK( 500. == Approx( resonance.capture().value ) );
      CHECK( 400. == Approx( resonance.fission().value ) );
      CHECK( 200. == Approx( resonance.competition().value ) );
    } // THEN
  } // GIVEN

  GIVEN( "a valid ResonanceTable where one channel has "
         "energy dependent widths" ) {

    ResonanceTable table( { { 1000. * electronVolt,
                              250. * electronVolt,
                              100. * rootElectronVolt,
                              500. * electronVolt,
                              500. * electronVolt,
                              300. * electronVolt },
                            { 5000. * electronVolt,
                              250. * electronVolt,
                              100. * rootElectronVolt,
                              500. * electronVolt,
                              1000. * electronVolt,
                              300. * electronVolt },
                            { 8000. * electronVolt,
                              250. * electronVolt,
                              100. * rootElectronVolt,
                              500. * electronVolt,
                              600. * electronVolt,
                              300. * electronVolt } },
                          { 2, 0, 4, 1 } );

    THEN( "Unresolved resonances can be retrieved at any energy" ) {

      Resonance resonance = table( 1000. * electronVolt );

      CHECK( 1000. == Approx( resonance.energy().value ) );
      CHECK( 250. == Approx( resonance.levelSpacing().value ) );
      CHECK( 100. == Approx( resonance.elastic().value ) );
      CHECK( 500. == Approx( resonance.capture().value ) );
      CHECK( 500. == Approx( resonance.fission().value ) );
      CHECK( 300. == Approx( resonance.competition().value ) );

      resonance = table( 3000. * electronVolt );

      CHECK( 3000. == Approx( resonance.energy().value ) );
      CHECK( 250. == Approx( resonance.levelSpacing().value ) );
      CHECK( 100. == Approx( resonance.elastic().value ) );
      CHECK( 500. == Approx( resonance.capture().value ) );
      CHECK( 750. == Approx( resonance.fission().value ) );
      CHECK( 300. == Approx( resonance.competition().value ) );

      resonance = table( 5000. * electronVolt );

      CHECK( 5000. == Approx( resonance.energy().value ) );
      CHECK( 250. == Approx( resonance.levelSpacing().value ) );
      CHECK( 100. == Approx( resonance.elastic().value ) );
      CHECK( 500. == Approx( resonance.capture().value ) );
      CHECK( 1000. == Approx( resonance.fission().value ) );
      CHECK( 300. == Approx( resonance.competition().value ) );

      resonance = table( 6500. * electronVolt );

      CHECK( 6500. == Approx( resonance.energy().value ) );
      CHECK( 250. == Approx( resonance.levelSpacing().value ) );
      CHECK( 100. == Approx( resonance.elastic().value ) );
      CHECK( 500. == Approx( resonance.capture().value ) );
      CHECK( 800. == Approx( resonance.fission().value ) );
      CHECK( 300. == Approx( resonance.competition().value ) );

      resonance = table( 8000. * electronVolt );

      CHECK( 8000. == Approx( resonance.energy().value ) );
      CHECK( 250. == Approx( resonance.levelSpacing().value ) );
      CHECK( 100. == Approx( resonance.elastic().value ) );
      CHECK( 500. == Approx( resonance.capture().value ) );
      CHECK( 600. == Approx( resonance.fission().value ) );
      CHECK( 300. == Approx( resonance.competition().value ) );
    } // THEN
  } // GIVEN

  GIVEN( "a valid ResonanceTable with full energy dependence" ) {

    ResonanceTable table( { { 1000. * electronVolt,
                              350. * electronVolt,
                              800. * rootElectronVolt,
                              500. * electronVolt,
                              100. * electronVolt,
                              200. * electronVolt },
                            { 5000. * electronVolt,
                              250. * electronVolt,
                              100. * rootElectronVolt,
                              1000. * electronVolt,
                              300. * electronVolt,
                              600. * electronVolt },
                            { 8000. * electronVolt,
                              550. * electronVolt,
                              200. * rootElectronVolt,
                              600. * electronVolt,
                              500. * electronVolt,
                              1000. * electronVolt } },
                            { 2, 0, 4, 1 } );

    THEN( "Unresolved resonances can be retrieved at any energy" ) {

      Resonance resonance = table( 1000. * electronVolt );

      CHECK( 1000. == Approx( resonance.energy().value ) );
      CHECK( 350. == Approx( resonance.levelSpacing().value ) );
      CHECK( 800. == Approx( resonance.elastic().value ) );
      CHECK( 500. == Approx( resonance.capture().value ) );
      CHECK( 100. == Approx( resonance.fission().value ) );
      CHECK( 200. == Approx( resonance.competition().value ) );

      resonance = table( 3000. * electronVolt );

      CHECK( 3000. == Approx( resonance.energy().value ) );
      CHECK( 300. == Approx( resonance.levelSpacing().value ) );
      CHECK( 450. == Approx( resonance.elastic().value ) );
      CHECK( 750. == Approx( resonance.capture().value ) );
      CHECK( 200. == Approx( resonance.fission().value ) );
      CHECK( 400. == Approx( resonance.competition().value ) );

      resonance = table( 5000. * electronVolt );

      CHECK( 5000. == Approx( resonance.energy().value ) );
      CHECK( 250. == Approx( resonance.levelSpacing().value ) );
      CHECK( 100. == Approx( resonance.elastic().value ) );
      CHECK( 1000. == Approx( resonance.capture().value ) );
      CHECK( 300. == Approx( resonance.fission().value ) );
      CHECK( 600. == Approx( resonance.competition().value ) );

      resonance = table( 6500. * electronVolt );

      CHECK( 6500. == Approx( resonance.energy().value ) );
      CHECK( 400. == Approx( resonance.levelSpacing().value ) );
      CHECK( 150. == Approx( resonance.elastic().value ) );
      CHECK( 800. == Approx( resonance.capture().value ) );
      CHECK( 400. == Approx( resonance.fission().value ) );
      CHECK( 800. == Approx( resonance.competition().value ) );

      resonance = table( 8000. * electronVolt );

      CHECK( 8000. == Approx( resonance.energy().value ) );
      CHECK( 550. == Approx( resonance.levelSpacing().value ) );
      CHECK( 200. == Approx( resonance.elastic().value ) );
      CHECK( 600. == Approx( resonance.capture().value ) );
      CHECK( 500. == Approx( resonance.fission().value ) );
      CHECK( 1000. == Approx( resonance.competition().value ) );
    } // THEN
  } // GIVEN

  GIVEN( "a valid ResonanceTable" ) {

    ResonanceTable table( { { 1000. * electronVolt,
                               250. * electronVolt,
                               100. * rootElectronVolt,
                               500. * electronVolt,
                               400. * electronVolt,
                               200. * electronVolt },
                            { 5000. * electronVolt,
                               250. * electronVolt,
                               100. * rootElectronVolt,
                               500. * electronVolt,
                               400. * electronVolt,
                               200. * electronVolt } },
                          { 2, 0, 4, 1 } );

    THEN( "an exception is thrown when going outside the energy range" ) {

      CHECK_THROWS( table( 999.9999 * electronVolt ) );
      CHECK_THROWS( table( 5000.0001 * electronVolt ) );
    } // THEN
  } // GIVEN
} // SCENARIO
