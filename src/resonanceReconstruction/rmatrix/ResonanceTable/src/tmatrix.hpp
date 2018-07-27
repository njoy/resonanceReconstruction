/**
 *  @brief Return T = ( 1 - RL )^-1 R matrix for the resonance table
 */
Matrix< std::complex< double > >
tmatrix( const Energy& energy,
         std::vector< std::complex< double > >& diagonalLMatrix ) const {

// BEGIN REALLY BAD FOR NOW - GET TESTING GOING
  // range with the R-matrices for each resonance
  auto rmatrices = this->resonances()
                     | ranges::view::transform(
                         [&] ( const auto& resonance )
                             { return resonance.rmatrix( energy ); } );

  // accumulate the rmatrix
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

  // calculate and return T = ( 1 - RL )^-1 R
  return
    ( Matrix< double >::Identity( size, size ) -
      rMatrix *
      Eigen::Map< Eigen::VectorXcd >( diagonalLMatrix.data(),
                                      size ).asDiagonal() ).inverse() *
    rMatrix;
// END REALLY BAD FOR NOW - GET TESTING GOING
}

