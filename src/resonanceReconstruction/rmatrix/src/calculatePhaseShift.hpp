/**
 *  @brief Default value for the phase shift
 *
 *  For all channel types except for neutron and charged particle channels, the 
 *  shift factor is 0.0.
 */
template < typename Type >
inline
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
inline
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
inline
double calculatePhaseShift< ChargedParticle >( const unsigned int,
                                               const double,
                                               const double ) {
  //! @todo develop this
  return 0.0;
}

