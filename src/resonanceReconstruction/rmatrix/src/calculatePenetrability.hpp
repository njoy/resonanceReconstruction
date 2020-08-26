/**
 *  @brief Default value for the penetrability
 *
 *  For all channel types except for neutron and charged particle channels, the
 *  penetrability is 1.0.
 */
template < typename Type >
double calculatePenetrability( const unsigned int, const double, const double ) {

  return 1.0;
}

/**
 *  @brief The penetrability for neutron channels
 *
 *  For a neutron channel, the penetrability is defined using the neutron
 *  hardsphere functions.
 *
 *  @param[in] l       the oribital angular momentum
 *  @param[in] ratio   the value of rho = ka
 */
template <>
double calculatePenetrability< Neutron >( const unsigned int l,
                                          const double ratio,
                                          const double  ) {

  return resonanceReconstruction::penetrationShift( l, ratio )[0];
}

/**
 *  @brief The penetrability for charged particle channels
 *
 *  For a charged particle channel, the penetrability is defined using the
 *  Coulomb wave functions.
 *
 *  @param[in] l       the oribital angular momentum
 *  @param[in] ratio   the value of rho = ka
 *  @param[in] eta     the Sommerfeld parameter
 */
template <>
double calculatePenetrability< ChargedParticle >( const unsigned int l,
                                                  const double ratio,
                                                  const double eta ) {


  std::complex< double > gf, dgf;
  coulombWaveFunctions( l, ratio, eta, gf, dgf);

  double F = gf.imag();
  double G = gf.real();
  return ( F == 0. ) and ( G == 0. ) ? 0. : ratio / ( F * F + G * G );
}
