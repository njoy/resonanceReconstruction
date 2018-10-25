/**
 *  @brief Default value for the shift factor
 *
 *  For all channel types except for neutron channels, the shift factor is 0.0.
 */
template < typename Type >
inline
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
inline
double calculateShiftFactor< Neutron >( const unsigned int l,
                                        const double ratio,
                                        const double ) {
  return resonanceReconstruction::penetrationShift( l, ratio )[1];
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
inline
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

