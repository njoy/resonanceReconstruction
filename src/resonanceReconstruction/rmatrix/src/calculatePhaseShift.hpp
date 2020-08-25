/**
 *  @brief Default value for the phase shift
 *
 *  For all channel types except for neutron and charged particle channels, the
 *  shift factor is 0.0.
 */
template < typename Type >
double calculatePhaseShift( const unsigned int, const double, const double ) {

  return 0.0;
}

/**
 *  @brief The phase shift for neutron channels
 *
 *  For a neutron channel, the phase shift is defined using the neutron
 *  hardsphere functions.
 *
 *  @param[in] l       the oribital angular momentum
 *  @param[in] ratio   the value of rho = ka
 */
template <>
double calculatePhaseShift< Neutron >( const unsigned int l,
                                       const double ratio,
                                       const double ) {

  return resonanceReconstruction::phaseShift( l, ratio );
}

/**
 *  @brief The phase shift for charged particle channels
 *
 *  For a charged particle channel, the phase shift is defined using the
 *  Coulomb wave functions.
 *
 *  @param[in] l       the oribital angular momentum
 *  @param[in] ratio   the value of rho = ka
 *  @param[in] eta     the Sommerfeld parameter
 */
template <>
double calculatePhaseShift< ChargedParticle >( const unsigned int l,
                                               const double ratio,
                                               const double eta ) {

  std::complex< double > gf, dgf;
  coulombWaveFunctions( l, ratio, eta, gf, dgf);

  double F = gf.imag();
  double G = gf.real();
  return ( F == 0. ) and ( G == 0. )
           ? 0.
           : std::acos( G / std::sqrt( F * F + G * G ) );
}
