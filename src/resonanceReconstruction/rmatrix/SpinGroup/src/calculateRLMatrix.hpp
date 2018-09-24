void calculateRLMatrix( const Energy& energy, ReichMoore ) {

// BEGIN REALLY BAD FOR NOW - GET TESTING GOING
  // range with the R-matrices for each resonance
  auto rmatrices = this->resonanceTable().resonances()
                     | ranges::view::transform(
                         [&] ( const auto& resonance )
                             { return resonance.rmatrix( energy ); } );

  // accumulate the rmatrix
  const unsigned int size = this->resonanceTable().numberChannels();
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
  this->matrix_ = Matrix< double >::Identity( size, size );
  this->matrix_ -= rMatrix *
             Eigen::Map< Eigen::VectorXcd >( this->diagonalLMatrix_.data(),
                                             size ).asDiagonal();
  this->matrix_ = this->matrix_.inverse();
  this->matrix_ *= rMatrix;
// END REALLY BAD FOR NOW - GET TESTING GOING
}

