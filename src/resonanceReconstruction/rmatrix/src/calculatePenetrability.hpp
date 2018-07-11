/**
 *  @brief Default value for the penetrability
 *
 *  For all channel types except for neutron channels, the penetrability is 1.0.
 */
template < typename Type >
inline
double calculatePenetrability( const unsigned int, const double ) {
  return 1.0;
}

/**
 *  @brief The penetrability for neutron channels
 *
 *  For a neutron channel, the penetrability is defined using the neutron 
 *  hardsphere functions.
 */
template <> 
inline
double calculatePenetrability< Neutron >( const unsigned int l,
                                          const double ratio ) {
  return resonanceReconstruction::penetrationShift( l, ratio )[0];
}

