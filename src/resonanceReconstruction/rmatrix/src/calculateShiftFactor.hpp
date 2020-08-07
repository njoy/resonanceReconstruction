/**
 *  @brief Default value for the shift factor
 *
 *  For all channel types except for neutron channels, the shift factor is 0.0.
 */
template < typename Type >
double calculateShiftFactor( const unsigned int, const double, const double ) {

  return 0.0;
}

/**
 *  @brief The penetrability for neutron channels
 *
 *  For a neutron channel, the shift factor is defined using the neutron
 *  hardsphere functions.
 */
template <>
double calculateShiftFactor< Neutron >( const unsigned int l,
                                        const double ratio,
                                        const double ) {

  double squared = ratio * ratio;
  switch ( l ) {

    case 0 : return 0.;
    case 1 : return - 1. / ( 1. + squared );
    case 2 : {

      constexpr std::array< double, 2 > numerator = {{ 18., 3. }};
      constexpr std::array< double, 3 > denominator = {{ 9., 3., 1. }};
      return -horner( std::rbegin( denominator ), std::rend( denominator ), squared )
             / horner( std::rbegin( denominator ), std::rend( denominator ), squared );
    }
    case 3 : {

      constexpr std::array< double, 3 > numerator = {{ 675., 90., 6. }};
      constexpr std::array< double, 4 > denominator = {{ 225., 45., 6., 1. }};
      return -horner( std::rbegin( denominator ), std::rend( denominator ), squared )
             / horner( std::rbegin( denominator ), std::rend( denominator ), squared );
    }
    case 4 : {

      constexpr std::array< double, 4 > numerator = {{ 44100., 4725., 270., 10. }};
      constexpr std::array< double, 5 > denominator = {{ 11025., 1575., 135., 10., 1. }};
      return -horner( std::rbegin( denominator ), std::rend( denominator ), squared )
             / horner( std::rbegin( denominator ), std::rend( denominator ), squared );
    }
    default : throw std::exception();
  }
}

/**
 *  @brief The shift factor for charged particle channels
 *
 *  For a charged particle channel, the shift factor is defined using the
 *  Coulomb wave functions.
 *
 *  @param[in] l       the oribital angular momentum
 *  @param[in] ratio   the value of rho = ka
 *  @param[in] eta     the Sommerfeld parameter
 */
template <>
double calculateShiftFactor< ChargedParticle >( const unsigned int l,
                                                const double ratio,
                                                const double eta ) {

  std::complex< double > gf, dgf;
  coulombWaveFunctions( l, ratio, eta, gf, dgf);

  double F = gf.imag();
  double G = gf.real();
  double dF = dgf.imag();
  double dG = dgf.real();
  return ratio / ( F * F + G * G ) * ( F * dF + G * dG );
}
