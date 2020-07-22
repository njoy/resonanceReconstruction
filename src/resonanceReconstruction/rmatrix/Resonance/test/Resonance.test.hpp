SCENARIO( "Resonance" ) {

  GIVEN( "valid data for a Resonance" ) {

    // single resonance data
    Energy energy = 6.823616e+4 * electronVolt;
    ReducedWidth eliminated = 3.933600e-1 * rootElectronVolt;
    std::vector< ReducedWidth > widths = { 2.179040e+2 * rootElectronVolt,
                                           1.000000e-5 * rootElectronVolt };

    THEN( "a Resonance can be constructed with an eliminated width" ) {

      Resonance resonance( energy, std::move( widths ), eliminated );

      CHECK( 6.823616e+4 == Approx( resonance.energy().value ) );

      CHECK( 3.933600e-1 == Approx( resonance.eliminatedWidth().value ) );

      CHECK( 2 == resonance.widths().size() );
      CHECK( 2.179040e+2 == Approx( resonance.widths()[0].value ) );
      CHECK( 1.000000e-5 == Approx( resonance.widths()[1].value ) );
    } // THEN

    THEN( "a Resonance can be constructed without an eliminated width" ) {

      Resonance resonance( energy, std::move( widths ) );

      CHECK( 6.823616e+4 == Approx( resonance.energy().value ) );

      CHECK( 0.0 == Approx( resonance.eliminatedWidth().value ) );

      CHECK( 2 == resonance.widths().size() );
      CHECK( 2.179040e+2 == Approx( resonance.widths()[0].value ) );
      CHECK( 1.000000e-5 == Approx( resonance.widths()[1].value ) );
    } // THEN
  } // GIVEN
} // SCENARIO
