/* D accounts for the statistical factors of "repeated" Jvalues for use in
 * computing the potential scattering.
 *
 * Its value can be expressed as:
 *
 *                  I + l + 1/2
 *  D = 2l + 1 -        Σ gj
 *               J = ||I - l| - 1/2|
 *
 * Noting gj is defined as
 *
 *  gj = (2J + 1) / (4I + 2)
 *
 * D can be interpretted as a function of I and l. Given l can assume a value
 * between 0 and 4 inclusively in ENDF and I is always a non-negative multiple
 * of 1/2, these values may be precomputed thereby avoid the evaluation of a
 * loop at runtime.
 *
 * The value of D for a given l value converges as I becomes large.
 * For values of I > 7/2,
 *
 * D(I,l) = l
 */
inline double D( const double I, const int l ){
  static constexpr const std::array< double, 40 > cache =
  {{ 0., 0., 0., 0., 0.,
     0., 0.75, 1.25, 1.75, 2.25 ,
     0., 1., 1.666666666666667, 2.333333333333333, 3.,
     0., 1., 1.875, 2.625, 3.375,
     0., 1., 2., 2.8, 3.6,
     0., 1., 2., 2.916666666666667, 3.75,
     0., 1., 2., 3., 3.857142857142858,
     0., 1., 2., 3., 3.9375 }};

  return ( I > 3.5 ) ? l : cache[ std::lround( 10. * I ) + l ];
}
