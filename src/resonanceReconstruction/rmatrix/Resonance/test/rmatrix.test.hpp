SCENARIO( "Resonance.rmatrix( energy )" ) {

  GIVEN( "Valid Resonances" ) {

    Resonance resonance1( 6.823616e+4 * electronVolt,
                         { 2.179040e+2 * rootElectronVolt,
                           1.000000e-5 * rootElectronVolt },
                         3.933600e-1 * rootElectronVolt );

    Resonance resonance2( 1.150980e+5 * electronVolt,
                          { 4.307780e+0 * rootElectronVolt,
                            0.0 * rootElectronVolt },
                          7.390000e-1 * rootElectronVolt );

    Resonance resonance3( 1.825230e+5 * electronVolt,
                         { 1.759740e+3 * rootElectronVolt,
                           4.000000e-1 * rootElectronVolt },
                         7.451500e-1 * rootElectronVolt );

    THEN( "an R-matrix contribution can be calculated for all of them" ) {

      auto rmatrix = resonance1.rmatrix( 1e+4 * electronVolt );

      CHECK( rmatrix.size() == 2 );
      CHECK( rmatrix.front().size() == 2 );

      CHECK( 8.153380E-01 == Approx( rmatrix[0][0].real() ) );
      CHECK( 2.166334E-06 == Approx( rmatrix[0][0].imag() ) );
      CHECK( 3.741730E-08 == Approx( rmatrix[0][1].real() ) );
      CHECK( 9.941688E-14 == Approx( rmatrix[0][1].imag() ) );
      CHECK( 3.741730E-08 == Approx( rmatrix[1][0].real() ) );
      CHECK( 9.941688E-14 == Approx( rmatrix[1][0].imag() ) );
      CHECK( 1.717146E-15 == Approx( rmatrix[1][1].real() ) );
      CHECK( 1.717146E-15 == Approx( rmatrix[1][1].imag() ) );

      rmatrix = resonance2.rmatrix( 1e+4 * electronVolt );

      CHECK( rmatrix.size() == 2 );
      CHECK( rmatrix.front().size() == 2 );

      CHECK( 1.765682E-04 == Approx( rmatrix[0][0].real() ) );
      CHECK( 9.175020E-10 == Approx( rmatrix[0][0].imag() ) );
      CHECK( 0.0 == Approx( rmatrix[0][1].real() ) );
      CHECK( 0.0 == Approx( rmatrix[0][1].imag() ) );
      CHECK( 0.0 == Approx( rmatrix[1][0].real() ) );
      CHECK( 0.0 == Approx( rmatrix[1][0].imag() ) );
      CHECK( 0.0 == Approx( rmatrix[1][1].real() ) );
      CHECK( 0.0 == Approx( rmatrix[1][1].imag() ) );

      rmatrix = resonance3.rmatrix( 1e+4 * electronVolt );

      CHECK( rmatrix.size() == 2 );
      CHECK( rmatrix.front().size() == 2 );

      CHECK( 1.794940E+01 == Approx( rmatrix[0][0].real() ) );
      CHECK( 5.776841E-05 == Approx( rmatrix[0][0].imag() ) );
      CHECK( 4.080013E-03 == Approx( rmatrix[0][1].real() ) );
      CHECK( 1.313112E-08 == Approx( rmatrix[0][1].imag() ) );
      CHECK( 4.080013E-03 == Approx( rmatrix[1][0].real() ) );
      CHECK( 1.313112E-08 == Approx( rmatrix[1][0].imag() ) );
      CHECK( 9.274126E-07 == Approx( rmatrix[1][1].real() ) );
      CHECK( 9.274126E-07 == Approx( rmatrix[1][1].imag() ) );
    }
  } // GIVEN
} // SCENARIO
