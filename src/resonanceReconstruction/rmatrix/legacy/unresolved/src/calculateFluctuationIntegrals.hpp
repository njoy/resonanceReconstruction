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

    {{ 1.1120413e-1, 2.3546798e-1, 2.8440987e-1, 2.2419127e-1 }},
    {{ 0.10967668e0, .030493789e0, 0.0042930874e0, 2.5827047e-4 }},
    {{ 4.9031965e-6, 1.4079206e-8, 0.033773418e0, 0.079932171e0 }},
    {{ 0.12835937e0, 0.17652616e0, 0.21347043e0, 0.21154965e0 }},
    {{ 0.13365186e0, 0.022630659e0, 1.6313638e-5, 2.745383e-31 }},
    {{ 3.3376214e-4, 0.018506108e0, 0.12309946e0, 0.29918923e0 }},
    {{ 0.33431475e0, 0.17766657e0, 0.042695894e0, 4.0760575e-3 }},
    {{ 1.1766115e-4, 5.0989546e-7, 1.7623788e-3, 0.021517749e0 }},
    {{ 0.080979849e0, 0.18797998e0, 0.30156335e0, 0.29616091e0 }},
    {{ 0.10775649e0, 2.5171914e-3, 8.9630388e-10, 0.e0 }}
  }};

  static constexpr std::array< std::array< double, 4 >, 10 > q = {{

    {{ 3.0013465e-3, 7.8592886e-2, 4.3282415e-1, 1.3345267e0 }},
    {{ 3.0481846e0, 5.8263198e0, 9.9452656e0, 1.5782128e1 }},
    {{ 23.996824e0, 36.216208e0, 1.3219203e-2, 7.2349624e-2 }},
    {{ 0.19089473e0, 0.39528842e0, 0.74083443e0, 1.3498293e0 }},
    {{ 2.5297983e0, 5.2384894e0, 13.821772e0, 75.647525e0 }},
    {{ 1.0004488e-3, 0.026197629e0, 0.14427472e0, 0.44484223e0 }},
    {{ 1.0160615e0, 1.9421066e0, 3.3150885e0, 5.2607092e0 }},
    {{ 7.9989414e0, 12.072069e0, 0.013219203e0, 0.072349624e0 }},
    {{ 0.19089473e0, 0.39528842e0, 0.74083443e0, 1.3498293e0 }},
    {{ 2.5297983e0, 5.2384894e0, 13.821772e0, 75.647525e0 }}
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
  else if ( not widths.hasCompetition() ) {

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
  else if ( not widths.hasFission() ) {

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
