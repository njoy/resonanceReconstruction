SCENARIO( "ResonanceTable" ) {

  GIVEN( "valid data for a ResonanceTable" ) {

    // single resonance data
    std::vector< Resonance > resonances =
      { Resonance( 6.823616e+4 * electronVolt,
                   3000. * electronVolt,
                   2.179040e+2 * rootElectronVolt,
                   1.000000e-5 * electronVolt,
                   20. * electronVolt,
                   5. * electronVolt ),
        Resonance( 1.150980e+5 * electronVolt,
                   2000. * electronVolt,
                   4.307780e+0 * rootElectronVolt,
                   0.0 * electronVolt,
                   4. * electronVolt,
                   6. * electronVolt ),
        Resonance( 1.825230e+5 * electronVolt,
                   1000. * electronVolt,
                   1.759740e+3 * rootElectronVolt,
                   4.000000e-1 * electronVolt,
                   1. * electronVolt,
                   3. * electronVolt ) };
    Degrees degrees = { 2, 0, 4, 1 };

    THEN( "a ResonanceTable can be constructed" ) {

      ResonanceTable table( std::move( resonances ),
                            std::move( degrees ) );

      CHECK( 2 == table.degreesOfFreedom().elastic );
      CHECK( 0 == table.degreesOfFreedom().capture );
      CHECK( 4 == table.degreesOfFreedom().fission );
      CHECK( 1 == table.degreesOfFreedom().competition );

      CHECK( 3 == table.numberResonances() );
      CHECK( 3 == table.resonances().size() );
      CHECK( 3 == table.energies().size() );
      CHECK( 6.823616e+4 == Approx( table.energies()[0].value ) );
      CHECK( 1.150980e+5 == Approx( table.energies()[1].value ) );
      CHECK( 1.825230e+5 == Approx( table.energies()[2].value ) );

      auto resonance = table.resonances()[0];
      CHECK( 6.823616e+4 == Approx( resonance.energy().value ) );
      CHECK( 3000. == Approx( resonance.levelSpacing().value ) );
      CHECK( 2.179040e+2 == Approx( resonance.elastic().value ) );
      CHECK( 1.000000e-5 == Approx( resonance.capture().value ) );
      CHECK( 20. == Approx( resonance.fission().value ) );
      CHECK( 5. == Approx( resonance.competition().value ) );

      resonance = table.resonances()[1];
      CHECK( 1.150980e+5 == Approx( resonance.energy().value ) );
      CHECK( 2000. == Approx( resonance.levelSpacing().value ) );
      CHECK( 4.307780e+0 == Approx( resonance.elastic().value ) );
      CHECK( 0.0 == Approx( resonance.capture().value ) );
      CHECK( 4. == Approx( resonance.fission().value ) );
      CHECK( 6. == Approx( resonance.competition().value ) );

      resonance = table.resonances()[2];
      CHECK( 1.825230e+5 == Approx( resonance.energy().value ) );
      CHECK( 1000. == Approx( resonance.levelSpacing().value ) );
      CHECK( 1.759740e+3 == Approx( resonance.elastic().value ) );
      CHECK( 4.000000e-1 == Approx( resonance.capture().value ) );
      CHECK( 1. == Approx( resonance.fission().value ) );
      CHECK( 3. == Approx( resonance.competition().value ) );
    } // THEN
  } // GIVEN

  GIVEN( "data for a ResonanceTable containing errors" ) {

    // wrong number of widths in a Resonance
    std::vector< Resonance > resonances =
      { Resonance( 6.823616e+4 * electronVolt,
                   3000. * electronVolt,
                   2.179040e+2 * rootElectronVolt,
                   1.000000e-5 * electronVolt,
                   20. * electronVolt,
                   5. * electronVolt ),
        Resonance( 1.150980e+5 * electronVolt,
                   2000. * electronVolt,
                   4.307780e+0 * rootElectronVolt,
                   0.0 * electronVolt,
                   4. * electronVolt,
                   6. * electronVolt ),
        Resonance( 1.825230e+5 * electronVolt,
                   1000. * electronVolt,
                   1.759740e+3 * rootElectronVolt,
                   4.000000e-1 * electronVolt,
                   1. * electronVolt,
                   3. * electronVolt ) };
    Degrees wrongValueDegrees = { 5, 0, 4, 1 };

    THEN( "an exception is thrown at construction when the value of degrees "
          "of freedom are not correct" ) {

      CHECK_THROWS( ResonanceTable( std::move( resonances ),
                                    std::move( wrongValueDegrees ) ) );
    } // THEN

    THEN( "an exception is thrown at construction when there are less than two "
          "resonances" ) {

      CHECK_THROWS( ResonanceTable( {}, { 1, 0, 1, 1 } ) );
      CHECK_THROWS( ResonanceTable(
                      { Resonance( 1e+3 * electronVolt,
                                   1e+1 * electronVolt,
                                   1e+2 * rootElectronVolt,
                                   1e-5 * electronVolt,
                                   20. * electronVolt,
                                   5. * electronVolt ) }, { 1, 0, 1, 1 } ) );
    } // THEN
  } // GIVEN
} // SCENARIO
