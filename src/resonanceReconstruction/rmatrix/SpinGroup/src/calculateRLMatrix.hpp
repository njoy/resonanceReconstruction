void calculateRLMatrix( const Energy& energy, ReichMoore ) {

  // range with the R-matrices for each resonance
  auto rmatrices = this->resonanceTable().resonances()
                     | ranges::view::transform(
                         [&] ( const auto& resonance )
                             { return resonance.rmatrix( energy ); } );

  // accumulate the rmatrix
  const unsigned int size = this->resonanceTable().numberChannels();
  this->matrix_ = Matrix< std::complex< double > >::Zero( size, size );
  for ( const auto& rmatrix : rmatrices ) {
    for ( unsigned int c = 0; c < size; ++c ) {
      for ( unsigned int cprime = 0; cprime < size; ++cprime ) {
        this->matrix_( c, cprime ) += rmatrix[c][cprime];
      }
    }
  }

  // zero out threshold reactions
  auto belowThreshold = this->belowThreshold( energy );
  for ( unsigned int c = 0; c < size; ++c ) {

    if ( belowThreshold[c] ) {

      this->matrix_.row(c).setZero();
      this->matrix_.col(c).setZero();
    }
  }

  // calculate and return R_L = ( 1 - RL )^-1 R
  this->rlmatrix_ = Matrix< double >::Identity( size, size );
  this->rlmatrix_ -= this->matrix_ *
             Eigen::Map< Eigen::VectorXcd >( this->diagonalLMatrix_.data(),
                                             size ).asDiagonal();
  this->rlmatrix_ = this->rlmatrix_.inverse();
  this->rlmatrix_ *= this->matrix_;
}
