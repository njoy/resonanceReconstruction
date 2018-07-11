/**
 *  @brief Default value for the phase shift
 *
 *  For all channel types except for neutron and charged particle channels, the 
 *  shift factor is 0.0.
 */
template < typename Type >
inline
double calculatePhaseShift( const unsigned int, const double ) {
  return 0.0;
}

/**
 *  @brief The phase shift for neutron channels
 *
 *  For a neutron channel, the phase shift is defined using the neutron 
 *  hardsphere functions.
 */
template <> 
inline
double calculatePhaseShift< Neutron >( const unsigned int l,
                                       const double ratio ) {
  return resonanceReconstruction::phaseShift( l, ratio );
}

/**
 *  @brief The phase shift for charged particle channels
 *
 *  For a charged particle channel, the phase shift is defined using the 
 *  Coulomb wave functions.
 */
template <> 
inline
double calculatePhaseShift< ChargedParticle >( const unsigned int,
                                               const double ) {
  //! @todo develop this
  return 0.0;
}

