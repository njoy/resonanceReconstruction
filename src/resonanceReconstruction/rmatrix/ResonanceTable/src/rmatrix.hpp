/**
 *  @brief Return the rmatrix for the resonance table
 */
Matrix< std::complex< double > > rmatrix( const Energy& energy ) const {

  // range with the R-matrices for each resonance
  auto rmatrices = this->resonances()
                     | ranges::view::transform(
                         [&] ( const auto& resonance )
                             { return resonance.rmatrix( energy ); } );

// BEGIN REALLY BAD - GET TESTING GOING
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
// END REALLY BAD - GET TESTING GOING
}