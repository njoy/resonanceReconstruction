template < typename Type >
inline
double calculateCoulombPhaseShift( const unsigned int, const double ) {
  return 0.0;
}

template <> 
inline
double calculateCoulombPhaseShift< ChargedParticle >( const unsigned int l,
                                                      const double eta ) {
  if ( l == 0 ) {
    return 0.0;
  }
  else {
    return eta; //! @todo implement Wc formula for charged particles
  }
}

