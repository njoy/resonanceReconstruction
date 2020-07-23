SCENARIO( "call operator" ) {

  GIVEN( "a valid UnresolvedResonanceTable for energy independent parameters" ) {

    UnresolvedResonanceTable table( { "1", "2" },
                                    { { 1000. * electronVolt,
                                        250. * electronVolt,
                                        { 100. * rootElectronVolt,
                                          500. * rootElectronVolt } },
                                      { 5000. * electronVolt,
                                        250. * electronVolt,
                                        { 100. * rootElectronVolt,
                                          500. * rootElectronVolt } } },
                                    { 4, 0 } );

    THEN( "Unresolved resonances can be retrieved at any energy" ) {

      UnresolvedResonance resonance = table( 1000. * electronVolt );

      CHECK( 1000. == Approx( resonance.energy().value ) );
      CHECK( 250. == Approx( resonance.levelSpacing().value ) );
      CHECK( 2 == resonance.widths().size() );
      CHECK( 100. == Approx( resonance.widths()[0].value ) );
      CHECK( 500. == Approx( resonance.widths()[1].value ) );

      resonance = table( 3000. * electronVolt );

      CHECK( 3000. == Approx( resonance.energy().value ) );
      CHECK( 250. == Approx( resonance.levelSpacing().value ) );
      CHECK( 2 == resonance.widths().size() );
      CHECK( 100. == Approx( resonance.widths()[0].value ) );
      CHECK( 500. == Approx( resonance.widths()[1].value ) );

      resonance = table( 5000. * electronVolt );

      CHECK( 5000. == Approx( resonance.energy().value ) );
      CHECK( 250. == Approx( resonance.levelSpacing().value ) );
      CHECK( 2 == resonance.widths().size() );
      CHECK( 100. == Approx( resonance.widths()[0].value ) );
      CHECK( 500. == Approx( resonance.widths()[1].value ) );
    } // THEN
  } // GIVEN

  GIVEN( "a valid UnresolvedResonanceTable where one channel has "
         "energy dependent widths" ) {

    UnresolvedResonanceTable table( { "1", "2" },
                                    { { 1000. * electronVolt,
                                        250. * electronVolt,
                                        { 100. * rootElectronVolt,
                                          500. * rootElectronVolt } },
                                      { 5000. * electronVolt,
                                        250. * electronVolt,
                                        { 100. * rootElectronVolt,
                                          1000. * rootElectronVolt } },
                                      { 8000. * electronVolt,
                                        250. * electronVolt,
                                        { 100. * rootElectronVolt,
                                          600. * rootElectronVolt } } },
                                    { 4, 0 } );

    THEN( "Unresolved resonances can be retrieved at any energy" ) {

      UnresolvedResonance resonance = table( 1000. * electronVolt );

      CHECK( 1000. == Approx( resonance.energy().value ) );
      CHECK( 250. == Approx( resonance.levelSpacing().value ) );
      CHECK( 2 == resonance.widths().size() );
      CHECK( 100. == Approx( resonance.widths()[0].value ) );
      CHECK( 500. == Approx( resonance.widths()[1].value ) );

      resonance = table( 3000. * electronVolt );

      CHECK( 3000. == Approx( resonance.energy().value ) );
      CHECK( 250. == Approx( resonance.levelSpacing().value ) );
      CHECK( 2 == resonance.widths().size() );
      CHECK( 100. == Approx( resonance.widths()[0].value ) );
      CHECK( 750. == Approx( resonance.widths()[1].value ) );

      resonance = table( 5000. * electronVolt );

      CHECK( 5000. == Approx( resonance.energy().value ) );
      CHECK( 250. == Approx( resonance.levelSpacing().value ) );
      CHECK( 2 == resonance.widths().size() );
      CHECK( 100. == Approx( resonance.widths()[0].value ) );
      CHECK( 1000. == Approx( resonance.widths()[1].value ) );

      resonance = table( 6500. * electronVolt );

      CHECK( 6500. == Approx( resonance.energy().value ) );
      CHECK( 250. == Approx( resonance.levelSpacing().value ) );
      CHECK( 2 == resonance.widths().size() );
      CHECK( 100. == Approx( resonance.widths()[0].value ) );
      CHECK( 800. == Approx( resonance.widths()[1].value ) );

      resonance = table( 8000. * electronVolt );

      CHECK( 8000. == Approx( resonance.energy().value ) );
      CHECK( 250. == Approx( resonance.levelSpacing().value ) );
      CHECK( 2 == resonance.widths().size() );
      CHECK( 100. == Approx( resonance.widths()[0].value ) );
      CHECK( 600. == Approx( resonance.widths()[1].value ) );
    } // THEN
  } // GIVEN

  GIVEN( "a valid UnresolvedResonanceTable with full energy dependence" ) {

    UnresolvedResonanceTable table( { "1", "2" },
                                    { { 1000. * electronVolt,
                                        350. * electronVolt,
                                        { 800. * rootElectronVolt,
                                          500. * rootElectronVolt } },
                                      { 5000. * electronVolt,
                                        250. * electronVolt,
                                        { 100. * rootElectronVolt,
                                          1000. * rootElectronVolt } },
                                      { 8000. * electronVolt,
                                        550. * electronVolt,
                                        { 200. * rootElectronVolt,
                                          600. * rootElectronVolt } } },
                                    { 4, 0 } );

    THEN( "Unresolved resonances can be retrieved at any energy" ) {

      UnresolvedResonance resonance = table( 1000. * electronVolt );

      CHECK( 1000. == Approx( resonance.energy().value ) );
      CHECK( 350. == Approx( resonance.levelSpacing().value ) );
      CHECK( 2 == resonance.widths().size() );
      CHECK( 800. == Approx( resonance.widths()[0].value ) );
      CHECK( 500. == Approx( resonance.widths()[1].value ) );

      resonance = table( 3000. * electronVolt );

      CHECK( 3000. == Approx( resonance.energy().value ) );
      CHECK( 300. == Approx( resonance.levelSpacing().value ) );
      CHECK( 2 == resonance.widths().size() );
      CHECK( 450. == Approx( resonance.widths()[0].value ) );
      CHECK( 750. == Approx( resonance.widths()[1].value ) );

      resonance = table( 5000. * electronVolt );

      CHECK( 5000. == Approx( resonance.energy().value ) );
      CHECK( 250. == Approx( resonance.levelSpacing().value ) );
      CHECK( 2 == resonance.widths().size() );
      CHECK( 100. == Approx( resonance.widths()[0].value ) );
      CHECK( 1000. == Approx( resonance.widths()[1].value ) );

      resonance = table( 6500. * electronVolt );

      CHECK( 6500. == Approx( resonance.energy().value ) );
      CHECK( 400. == Approx( resonance.levelSpacing().value ) );
      CHECK( 2 == resonance.widths().size() );
      CHECK( 150. == Approx( resonance.widths()[0].value ) );
      CHECK( 800. == Approx( resonance.widths()[1].value ) );

      resonance = table( 8000. * electronVolt );

      CHECK( 8000. == Approx( resonance.energy().value ) );
      CHECK( 550. == Approx( resonance.levelSpacing().value ) );
      CHECK( 2 == resonance.widths().size() );
      CHECK( 200. == Approx( resonance.widths()[0].value ) );
      CHECK( 600. == Approx( resonance.widths()[1].value ) );
    } // THEN
  } // GIVEN

  GIVEN( "a valid UnresolvedResonanceTable" ) {

    UnresolvedResonanceTable table( { "1", "2" },
                                    { { 1000. * electronVolt,
                                        250. * electronVolt,
                                        { 100. * rootElectronVolt,
                                          500. * rootElectronVolt } },
                                      { 5000. * electronVolt,
                                        250. * electronVolt,
                                        { 100. * rootElectronVolt,
                                          500. * rootElectronVolt } } },
                                    { 4, 0 } );

    THEN( "an exception is thrown when going outside the energy range" ) {

      CHECK_THROWS( table( 999.9999 * electronVolt ) );
      CHECK_THROWS( table( 5000.0001 * electronVolt ) );
    } // THEN
  } // GIVEN
} // SCENARIO
