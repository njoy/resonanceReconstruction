template < typename Type >
inline
double calculatePenetrability( const unsigned int, const double ) {
  return 1.0;
}

template <> 
inline
double calculatePenetrability< Neutron >( const unsigned int l,
                                          const double ratio ) {
  return resonanceReconstruction::penetrationShift( l, ratio )[0];
}

