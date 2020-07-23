SCENARIO( "calculateFluctuationIntegrals" ) {

  GIVEN( "valid values for the widths" ) {

    THEN( "the appropriate fluctuation integrals are calculated "
          "when all widths are non-zero" ) {

      Degrees degrees( 1, 0, 2, 3 );
      Widths widths( 4.75405e-2 * electronVolts, 4.07e-2 * electronVolts,
                     2.842 * electronVolts, 0.402 * electronVolts );
      FluctuationIntegrals integrals =
      calculateFluctuationIntegrals( widths, degrees );

      CHECK( 1.5667666623387133 == Approx( integrals.elastic.value ) );
      CHECK( 0.56993816408468989 == Approx( integrals.capture.value ) );
      CHECK( 0.25126381048572832 == Approx( integrals.fission.value ) );
      CHECK( 0.4680023731 == Approx( integrals.competition.value ) );

      // REMARK: NJOY2016 does not calculate the fluctuation integral for
      //         the competitive reaction
    } // THEN

    THEN( "the appropriate fluctuation integrals are calculated "
          "when all widths are non-zero except for fission" ) {

      Degrees degrees( 1, 0, 0, 3 );
      Widths widths( 4.75405e-2 * electronVolts, 4.07e-2 * electronVolts,
                     0.0 * electronVolts, 2.842 * electronVolts );
      FluctuationIntegrals integrals =
      calculateFluctuationIntegrals( widths, degrees );

      CHECK( 1.7591280230875452 == Approx( integrals.elastic.value ) );
      CHECK( 0.65378100220430613 == Approx( integrals.capture.value ) );
      CHECK( 0. == Approx( integrals.fission.value ) );
      CHECK( 0.3130757611 == Approx( integrals.competition.value ) );

      // REMARK: NJOY2016 does not calculate the fluctuation integral for
      //         the competitive reaction
    } // THEN

    THEN( "the appropriate fluctuation integrals are calculated "
          "when all widths are non-zero except for competition" ) {

      // values tested against NJOY2016 using ENDF/B-VIII.0 Pu239
      // at 2500 eV

      Degrees degrees( 1, 0, 2, 0 );
      Widths widths( 4.75405e-2 * electronVolts, 4.07e-2 * electronVolts,
                     2.842 * electronVolts, 0.0 * electronVolts );
      FluctuationIntegrals integrals =
      calculateFluctuationIntegrals( widths, degrees );

      CHECK( 2.2806310870217965 == Approx( integrals.elastic.value ) );
      CHECK( 0.89457720997275580 == Approx( integrals.capture.value ) );
      CHECK( 0.30087206891891111 == Approx( integrals.fission.value ) );
      CHECK( 0. == Approx( integrals.competition.value ) );
    } // THEN

    THEN( "the appropriate fluctuation integrals are calculated "
          "when only elastic and capture widths are non-zero" ) {

      Degrees degrees( 1, 0, 0, 0 );
      Widths widths( 4.75405e-2 * electronVolts, 4.07e-2 * electronVolts,
                     0.0 * electronVolts, 0.0 * electronVolts );
      FluctuationIntegrals integrals =
      calculateFluctuationIntegrals( widths, degrees );

      CHECK( 14.395389203673844 == Approx( integrals.elastic.value ) );
      CHECK( 7.7551848262740855 == Approx( integrals.capture.value ) );
      CHECK( 0. == Approx( integrals.fission.value ) );
      CHECK( 0. == Approx( integrals.competition.value ) );
    } // THEN
  } // GIVEN
} // SCENARIO
