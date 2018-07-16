/**
 *  @brief Return the rmatrix for the resonance table
 */
Matrix< std::complex< double > > rmatrix( const Energy& energy ) const {

  // initialise the rmatrix
  const unsigned int size = this->numberChannels();
  Matrix< std::complex< double > > rMatrix( size, size );

  // range with the R-matrices for each resonance
  auto rmatrices = this->resonances()
                     | ranges::view::transform(
                         [&] ( const auto& resonance )
                             { return resonance.rmatrix( energy ); } );

// BEGIN REALLY BAD - GET TESTING GOING
  for ( const auto& rmatrix : rmatrices ) {
    for ( unsigned int c = 0; c < size; ++c ) {
      for ( unsigned int cprime = 0; cprime < size; ++cprime ) {
        rMatrix( c, cprime ) += rmatrix[c][cprime];
      }
    }
  }
// END REALLY BAD - GET TESTING GOING

/*  // some labdas to sum matrices
  auto sumRow =
    [] ( auto&& leftRow, const auto&& rightRow ) {
      ranges::copy( ranges::view::zip_with( hana::plus, leftRow, rightRow ),
                    leftRow );
      return leftRow;
    };
  auto sumMatrix =
    [=] ( auto&& leftMatrix, const auto&& rightMatrix ) {
      ranges::for_each(
        ranges::view::zip_with( sumRow, leftMatrix, rightMatrix ),
        [] ( auto&&... ) { return; } );
      return leftMatrix;
    };

  // range with the R-matrices for each resonance
  auto rmatrices = this->resonances()
                     | ranges::view::transform(
                         [&] ( const auto& resonance )
                             { return resonance.rmatrix( energy ); } );

  // accumulate the R-matrix
  auto matrixView =
    ranges::view::cartesian_product( ranges::view::indices( 0, size ),
                                     ranges::view::indices( 0, size ) )
      | ranges::view::transform(
          [&] ( const auto& indices ) -> decltype( auto )
              { return rMatrix( std::get< 0 >( indices ), 
                                std::get< 1 >( indices ) ); } )
      | ranges::view::chunk( size );
  ranges::accumulate( rmatrices, matrixView, sumMatrix );*/

  return rMatrix;
}

