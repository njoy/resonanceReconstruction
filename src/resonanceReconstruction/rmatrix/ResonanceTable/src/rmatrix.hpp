/**
 *  @brief Return the rmatrix for the resonance table
 *
 *  @param[in] energy   the energy at which the rmatrix must be evaluated
 */
Matrix< std::complex< double > > rmatrix( const Energy& energy ) const {

  //! @todo The rmatrix function on the resonance assumes Reich-Moore,
  //!       full rmatrix will require a templated function.

  // range with the R-matrices for each resonance
  auto rmatrices = this->resonances()
                     | ranges::view::transform(
                         [&] ( const auto& resonance )
                             { return resonance.rmatrix( energy ); } );

  const unsigned int size = this->numberChannels();
  Matrix< std::complex< double > > rMatrix =
      Matrix< std::complex< double > >::Zero( size, size );
  for ( const auto& rmatrix : rmatrices ) {
    for ( unsigned int c = 0; c < size; ++c ) {
      for ( unsigned int cprime = 0; cprime < size; ++cprime ) {
        rMatrix( c, cprime ) += rmatrix[c][cprime];
      }
    }
  }
  return rMatrix;
}
