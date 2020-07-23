/**
 *  @brief Calculate the fluctuation integrals for the legacy unresolved
 *         resonance
 *
 *  @param widths    the full widths for which to calculate the integrals
 *  @param degrees   the degrees of freedom for each width
 *
 *  @return the fluctuation integral data
 */
FluctuationIntegrals calculateFluctuationIntegrals( const Widths& widths,
                                                    const Degrees& degrees ) {

  static constexpr std::array< std::array< double, 4 >, 10 > w = {{

    {{ 1.1120413E-01, 3.3773418E-02, 3.3376214E-04, 1.7623788E-03 }},
    {{ 2.3546798E-01, 7.9932171E-02, 1.8506108E-02, 2.1517749E-02 }},
    {{ 2.8440987E-01, 1.2835937E-01, 1.2309946E-01, 8.0979849E-02 }},
    {{ 2.2419127E-01, 1.7652616E-01, 2.9918923E-01, 1.8797998E-01 }},
    {{ 1.0967668E-01, 2.1347043E-01, 3.3431475E-01, 3.0156335E-01 }},
    {{ 3.0493789E-02, 2.1154965E-01, 1.7766657E-01, 2.9616091E-01 }},
    {{ 4.2930874E-03, 1.3365186E-01, 4.2695894E-02, 1.0775649E-01 }},
    {{ 2.5827047E-04, 2.2630659E-02, 4.0760575E-03, 2.5171914E-03 }},
    {{ 4.9031965E-06, 1.6313638E-05, 1.1766115E-04, 8.9630388E-10 }},
    {{ 1.4079206E-08, 2.7453830E-31, 5.0989546E-07, 0.0000000E+00 }}
  }};

  static constexpr std::array< std::array< double, 4 >, 10 > q = {{

    {{ 3.0013465E-03, 1.3219203E-02, 1.0004488E-03, 1.3219203E-02 }},
    {{ 7.8592886E-02, 7.2349624E-02, 2.6197629E-02, 7.2349624E-02 }},
    {{ 4.3282415E-01, 1.9089473E-01, 1.4427472E-01, 1.9089473E-01 }},
    {{ 1.3345267E+00, 3.9528842E-01, 4.4484223E-01, 3.9528842E-01 }},
    {{ 3.0481846E+00, 7.4083443E-01, 1.0160615E+00, 7.4083443E-01 }},
    {{ 5.8263198E+00, 1.3498293E+00, 1.9421066E+00, 1.3498293E+00 }},
    {{ 9.9452656E+00, 2.5297983E+00, 3.3150885E+00, 2.5297983E+00 }},
    {{ 1.5782128E+01, 5.2384894E+00, 5.2607092E+00, 5.2384894E+00 }},
    {{ 2.3996824E+01, 1.3821772E+01, 7.9989414E+00, 1.3821772E+01 }},
    {{ 3.6216208E+01, 7.5647525E+01, 1.2072069E+01, 7.5647525E+01 }}
  }};

  // initialise result
  FluctuationIntegrals result( 0. / electronVolt, 0. / electronVolt,
                               0. / electronVolt, 0. / electronVolt );

  // the degrees of freedom for each reaction
  unsigned int mu = degrees.elastic - 1;
  unsigned int nu = degrees.hasFission() ? degrees.fission - 1 : 0;
  unsigned int lambda = degrees.hasCompetition() ? degrees.competition - 1 : 0;

  //! @todo we don't really need the fluctuation integral for the competition

  // case 1: fission and competition are zero
  if ( ( not widths.hasFission() ) and ( not widths.hasCompetition() ) ) {

    for ( unsigned int i = 0; i < 10; ++i ) {

      double factor = w[i][mu];
      Width denominator = ( q[i][mu] * widths.elastic + widths.capture );
      result += FluctuationIntegrals(
                    q[i][mu] * q[i][mu] * factor / denominator,
                    q[i][mu] * factor / denominator,
                    0. / electronVolt, 0. / electronVolt );
    }
  }
  // case 2: only competition is zero
  else if ( ( widths.hasFission() ) and ( not widths.hasCompetition() ) ) {

    for ( unsigned int i = 0; i < 10; ++i ) {

      for ( unsigned int j = 0; j < 10; ++j ) {

        double factor = w[i][mu] * w[j][nu];
        Width denominator = ( q[i][mu] * widths.elastic + widths.capture
                            + q[j][nu] * widths.fission );
        result += FluctuationIntegrals(
                      q[i][mu] * q[i][mu] * factor / denominator,
                      q[i][mu] * factor / denominator,
                      q[i][mu] * q[j][nu] * factor / denominator,
                      0. / electronVolt );
      }
    }
  }
  // case 3: only fission is zero
  else if ( ( not widths.hasFission() ) and ( widths.hasCompetition() ) ) {

    for ( unsigned int i = 0; i < 10; ++i ) {

      for ( unsigned int k = 0; k < 10; ++k ) {

        double factor = w[i][mu] * w[k][lambda];
        Width denominator = ( q[i][mu] * widths.elastic + widths.capture
                            + q[k][lambda] * widths.competition );
        result += FluctuationIntegrals(
                      q[i][mu] * q[i][mu] * factor / denominator,
                      q[i][mu] * factor / denominator,
                      0. / electronVolt,
                      q[i][mu] * q[k][lambda] * factor / denominator );
      }
    }
  }
  // case 4: all non-zero
  else {

    for ( unsigned int i = 0; i < 10; ++i ) {

      for ( unsigned int j = 0; j < 10; ++j ) {

        for ( unsigned int k = 0; k < 10; ++k ) {

          double factor = w[i][mu] * w[j][nu] * w[k][lambda];
          Width denominator = ( q[i][mu] * widths.elastic + widths.capture
                              + q[j][nu] * widths.fission
                              + q[k][lambda] * widths.competition );
          result += FluctuationIntegrals(
                        q[i][mu] * q[i][mu] * factor / denominator,
                        q[i][mu] * factor / denominator,
                        q[i][mu] * q[j][nu] * factor / denominator,
                        q[i][mu] * q[k][lambda] * factor / denominator );
        }
      }
    }
  }

  return result;
}
