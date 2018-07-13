/**
 *  @class
 *  @brief A spin group corresponding to a J,pi value
 */
class SpinGroup {

  /* fields */
  std::vector< Channel > channels_;
  ResonanceTable parameters_;

public:

  /* constructor */

  auto resonanceTable() const { return this->parameters_; }

  auto generateCollisionMatrix( const Energy&  ) const {

    // this code calculates the R-matrix using one eliminated channel
/*
    // initialise the rmatrix
    const int size = this->resonanceTable().numberChannels();
    Matrix< std::complex< double > > rMatrix( size, size );

    // some labdas to sum matrices
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
    auto rmatrices = this->resonanceTable().resonances()
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
  }
};
