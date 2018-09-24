/**
 *  @brief L matrix diagonal for the ShiftFactor boundary condition
 *
 *  The penetrability P and shift factor S are potentially energy dependent
 *  functions (for neutron and charged particle channels) while the boundary
 *  condition is a real value independent of energy.
 *
 *  @param[in] energy            the energy value
 *  @param[in] penetrabilities   the channel penetrabilities
 *
 *  @return The complex number S - B + iP
 */
template < typename Range >
auto calculateLDiagonal( const Energy& energy,
                             Range penetrabilities,
                             Constant ) const {

  return ranges::view::zip_with(
             [] ( double S, double B, double P )
                { return std::complex< double >( S - B, P ); },
             this->shiftFactors( energy ),
             this->boundaryConditions(),
             penetrabilities );
}

/**
 *  @brief L matrix diagonal for the ShiftFactor boundary condition
 *
 *  The R-matrix analysis code SAMMY uses a different approach to the boundary
 *  condition. SAMMY effectively eliminates the real part of the L matrix 
 *  diagonal by setting the boundary condition B to be equal to the shift 
 *  factor S, making the bondary condition potentially an energy dependent
 *  quantity.
 *
 *  Using the SAMMY boundary condition will effectively ignore the value of the
 *  boundary condition of each channel.
 *
 *  @param[in] energy            the energy value
 *  @param[in] penetrabilities   the channel penetrabilities
 *
 *  @return The complex number iP
 */
template < typename Range >
auto calculateLDiagonal( const Energy&,
                             Range penetrabilities,
                             ShiftFactor ) const {

  return penetrabilities
             | ranges::view::transform(
                   [] ( double P )
                      { return std::complex< double >( 0.0, P ); } );
}
