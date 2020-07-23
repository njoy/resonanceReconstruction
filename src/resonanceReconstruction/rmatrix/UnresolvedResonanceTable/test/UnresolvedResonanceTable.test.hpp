SCENARIO( "UnresolvedResonanceTable" ) {

  GIVEN( "valid data for a UnresolvedResonanceTable" ) {

    // single resonance data
    std::vector< ChannelID > channels = { "1", "2" };
    std::vector< UnresolvedResonance > resonances =
      { UnresolvedResonance( 6.823616e+4 * electronVolt,
                             3000. * electronVolt,
                             { 2.179040e+2 * rootElectronVolt,
                               1.000000e-5 * rootElectronVolt } ),
        UnresolvedResonance( 1.150980e+5 * electronVolt,
                             2000. * electronVolt,
                             { 4.307780e+0 * rootElectronVolt,
                               0.0 * rootElectronVolt } ),
        UnresolvedResonance( 1.825230e+5 * electronVolt,
                             1000. * electronVolt,
                             { 1.759740e+3 * rootElectronVolt,
                               4.000000e-1 * rootElectronVolt } ) };
    std::vector< unsigned int > degrees = { 4, 1 };

    THEN( "a ResonanceTable can be constructed" ) {

      UnresolvedResonanceTable table( std::move( channels ),
                                      std::move( resonances ),
                                      std::move( degrees ) );

      CHECK( 2 == table.numberChannels() );
      CHECK( 2 == table.channels().size() );
      CHECK( "1" == table.channels()[0] );
      CHECK( "2" == table.channels()[1] );

      CHECK( 2 == table.degreesOfFreedom().size() );
      CHECK( 4 == table.degreesOfFreedom()[0] );
      CHECK( 1 == table.degreesOfFreedom()[1] );

      CHECK( 3 == table.numberResonances() );
      CHECK( 3 == table.resonances().size() );
      CHECK( 3 == table.energies().size() );
      CHECK( 6.823616e+4 == Approx( table.energies()[0].value ) );
      CHECK( 1.150980e+5 == Approx( table.energies()[1].value ) );
      CHECK( 1.825230e+5 == Approx( table.energies()[2].value ) );

      auto resonance = table.resonances()[0];
      CHECK( 6.823616e+4 == Approx( resonance.energy().value ) );
      CHECK( 3000. == Approx( resonance.levelSpacing().value ) );
      CHECK( 2 == resonance.widths().size() );
      CHECK( 2.179040e+2 == Approx( resonance.widths()[0].value ) );
      CHECK( 1.000000e-5 == Approx( resonance.widths()[1].value ) );

      resonance = table.resonances()[1];
      CHECK( 1.150980e+5 == Approx( resonance.energy().value ) );
      CHECK( 2000. == Approx( resonance.levelSpacing().value ) );
      CHECK( 2 == resonance.widths().size() );
      CHECK( 4.307780e+0 == Approx( resonance.widths()[0].value ) );
      CHECK( 0.0 == Approx( resonance.widths()[1].value ) );

      resonance = table.resonances()[2];
      CHECK( 1.825230e+5 == Approx( resonance.energy().value ) );
      CHECK( 1000. == Approx( resonance.levelSpacing().value ) );
      CHECK( 2 == resonance.widths().size() );
      CHECK( 1.759740e+3 == Approx( resonance.widths()[0].value ) );
      CHECK( 4.000000e-1 == Approx( resonance.widths()[1].value ) );
    } // THEN
  } // GIVEN

  GIVEN( "data for a ResonanceTable containing errors" ) {

    // wrong number of widths in a Resonance
    std::vector< ChannelID > channels = { "1", "2" };
    std::vector< UnresolvedResonance > resonances =
      { UnresolvedResonance( 6.823616e+4 * electronVolt,
                             3000. * electronVolt,
                             { 2.179040e+2 * rootElectronVolt,
                               1.000000e-5 * rootElectronVolt } ),
        UnresolvedResonance( 1.150980e+5 * electronVolt,
                             2000. * electronVolt,
                             { 4.307780e+0 * rootElectronVolt,
                               0.0 * rootElectronVolt } ),
        UnresolvedResonance( 1.825230e+5 * electronVolt,
                             1000. * electronVolt,
                             { 1.759740e+3 * rootElectronVolt,
                               4.000000e-1 * rootElectronVolt } ) };
    std::vector< UnresolvedResonance > wrongResonances =
      { UnresolvedResonance(
                   6.823616e+4 * electronVolt,
                   3000. * electronVolt,
                   { 2.179040e+2 * rootElectronVolt,
                     1.000000e-5 * rootElectronVolt } ),
        UnresolvedResonance(
                   1.150980e+5 * electronVolt,
                   2000. * electronVolt,
                   { 4.307780e+0 * rootElectronVolt } ),
        UnresolvedResonance(
                   1.825230e+5 * electronVolt,
                   1000. * electronVolt,
                   { 1.759740e+3 * rootElectronVolt,
                     4.000000e-1 * rootElectronVolt } ) };
    std::vector< unsigned int > degrees = { 4, 1 };
    std::vector< unsigned int > wrongDegrees = { 4 };
    std::vector< unsigned int > wrongValueDegrees = { 5, 1 };

    THEN( "an exception is thrown at construction when a resonance does not "
          "contain the right number of widths" ) {

      auto copyChannels = channels;
      auto copyDegrees = degrees;
      CHECK_THROWS(
          UnresolvedResonanceTable( std::move( copyChannels ),
                                    std::move( wrongResonances ),
                                    std::move( copyDegrees ) ) );
    } // THEN

    THEN( "an exception is thrown at construction when the number of degrees "
          "of freedom are not correct" ) {

      auto copyChannels = channels;
      auto copyResonances = resonances;
      CHECK_THROWS(
          UnresolvedResonanceTable( std::move( copyChannels ),
                                    std::move( copyResonances ),
                                    std::move( wrongDegrees ) ) );
    } // THEN

    THEN( "an exception is thrown at construction when the value of degrees "
          "of freedom are not correct" ) {

      CHECK_THROWS(
          UnresolvedResonanceTable( std::move( channels ),
                                    std::move( resonances ),
                                    std::move( wrongValueDegrees ) ) );
    } // THEN
  } // GIVEN
} // SCENARIO
