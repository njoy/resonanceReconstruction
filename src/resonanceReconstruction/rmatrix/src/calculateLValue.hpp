/**
 *  @brief Default L matrix diagonal calculated as S - B + iP
 *
 *  The penetrability P and shift factor S are potentially energy dependent
 *  functions (for neutron and charged particle channels) while the boundary
 *  condition is a real value independent of energy.
 *
 *  @param[in] S : the shift factor
 *  @param[in] B : the boundary condition
 *  @param[in] P : the penetrability
 *
 *  @return The complex number S - B + iP
 */
template < typename BoundaryOption >
inline
std::complex< double >
calculateLValue( const double S, const double B, const double P ) {

  return std::complex< double >( S - B, P );
}

/**
 *  @brief L matrix diagonal for the Sammy boundary condition B = S
 *
 *  The R-matrix analysis code Sammy uses a different approach to the boundary
 *  condition. Sammy effectively eliminates the real part of the L matrix 
 *  diagonal by setting the boundary condition B to be equal to the shift 
 *  factor S, making the bondary condition potentially an energy dependent
 *  quantity.
 *
 *  @param[in] P : the penetrability
 *
 *  @return The complex number iP
 */
template <>
inline
std::complex< double >
calculateLValue< Sammy >( const double, const double, const double P ) {

  return std::complex< double >( 0.0, P );
}
