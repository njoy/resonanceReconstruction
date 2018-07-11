/**
 *  @brief Default value for the shift factor
 *
 *  For all channel types except for neutron channels, the shift factor is 0.0.
 */
template < typename Type >
inline
double calculateShiftFactor( const unsigned int, const double ) {
  return 0.0;
}

/**
 *  @brief The penetrability for neutron channels
 *
 *  For a neutron channel, the penetrability is defined using the neutron 
 *  hardsphere functions.
 */
template <> 
inline
double calculateShiftFactor< Neutron >( const unsigned int l,
                                        const double ratio ) {
  return resonanceReconstruction::penetrationShift( l, ratio )[1];
}

