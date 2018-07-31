/**
 *  @brief Return T = ( 1 - RL )^-1 R matrix for the resonance table
 *
 *  @param[in] energy            the energy for the T matrix evaluation
 *  @param[in] diagonalLMatrix   the diagonal of the L matrix
 *  @param[in,out] tMatrix       the initial value of the T matrix (identity
 *                               matrix for incident neutrons), stores the final
 *                               T matrix
 */
void tmatrix( const Energy& energy,
              std::vector< std::complex< double > >& diagonalLMatrix,
              Matrix< std::complex< double > >& tMatrix ) const {

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
  tMatrix -= rMatrix *
             Eigen::Map< Eigen::VectorXcd >( diagonalLMatrix.data(),
                                      size ).asDiagonal();
  tMatrix = tMatrix.inverse();
  tMatrix *= rMatrix;
// END REALLY BAD FOR NOW - GET TESTING GOING
}

