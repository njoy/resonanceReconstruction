SCENARIO( "calculateFluctuationIntegrals" ) {

  GIVEN( "valid values for the widths" ) {

    THEN( "the appropriate fluctuation integrals are calculated "
          "when all widths are non-zero" ) {

      Degrees degrees( 1, 1, 1, 1 );
      Widths widths( 1.0 * electronVolts, 2.0 * electronVolts,
                     3.0 * electronVolts, 4.0 * electronVolts );
      FluctuationIntegrals integrals =
      calculateFluctuationIntegrals( widths, degrees );

      CHECK( 1. / electronVolts == integrals.elastic );
      CHECK( 2. / electronVolts == integrals.capture );
      CHECK( 3. / electronVolts == integrals.fission );
      CHECK( 4. / electronVolts == integrals.competition );
    } // THEN

    THEN( "the appropriate fluctuation integrals are calculated "
          "when all widths are non-zero except for fission" ) {

      Degrees degrees( 1, 1, 0, 1 );
      Widths widths( 1.0 * electronVolts, 2.0 * electronVolts,
                     0.0 * electronVolts, 4.0 * electronVolts );
      FluctuationIntegrals integrals =
      calculateFluctuationIntegrals( widths, degrees );

      CHECK( 1. / electronVolts == integrals.elastic );
      CHECK( 2. / electronVolts == integrals.capture );
      CHECK( 0. / electronVolts == integrals.fission );
      CHECK( 4. / electronVolts == integrals.competition );
    } // THEN

    THEN( "the appropriate fluctuation integrals are calculated "
          "when all widths are non-zero except for competition" ) {

      Degrees degrees( 1, 1, 1, 0 );
      Widths widths( 1.0 * electronVolts, 2.0 * electronVolts,
                     3.0 * electronVolts, 0.0 * electronVolts );
      FluctuationIntegrals integrals =
      calculateFluctuationIntegrals( widths, degrees );

      CHECK( 1. / electronVolts == integrals.elastic );
      CHECK( 2. / electronVolts == integrals.capture );
      CHECK( 3. / electronVolts == integrals.fission );
      CHECK( 0. / electronVolts == integrals.competition );
    } // THEN

    THEN( "the appropriate fluctuation integrals are calculated "
          "when only elastic and capture widths are non-zero" ) {

      Degrees degrees( 1, 1, 0, 0 );
      Widths widths( 1.0 * electronVolts, 2.0 * electronVolts,
                     0.0 * electronVolts, 0.0 * electronVolts );
      FluctuationIntegrals integrals =
      calculateFluctuationIntegrals( widths, degrees );

      CHECK( 1. / electronVolts == integrals.elastic );
      CHECK( 2. / electronVolts == integrals.capture );
      CHECK( 0. / electronVolts == integrals.fission );
      CHECK( 0. / electronVolts == integrals.competition );
    } // THEN
  } // GIVEN
} // SCENARIO
