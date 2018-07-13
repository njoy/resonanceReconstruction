/**
 *  @brief Generate the upper triangular matrix with reduced width products
 *
 *  An element of this matrix corresponds with the product of the reduced
 *  widths of two channels.
 * 
 *
 *  The matrix returned from this function is an upper triangular matrix only
 *  as the matrix is symmetrical. To obtain the full matrix, symmetry must
 *  be applied.
 *
 *  @param[in] widths       the reduced widths (in sqrt(eV))
 */
static Matrix< double >
generateProducts( std::vector< ReducedWidth >&& widths ) {

  const int size = widths.size();
  Matrix< double > result( size, size );

  auto assign =
    [&] ( const int i, const int j )
        { result( i, j ) = widths[i] * widths[j] / electronVolt; };

  ranges::for_each(
      ranges::view::indices( 0, size ),
      [&] ( const int c )
          { ranges::for_each(
                ranges::view::indices( c, size ),
                [&] ( const int cprime )
                    { assign( c, cprime ); } ); } );
  //
  // auto squareWidths = ranges::view::cartesian_product( widths, widths )
  //                     | ranges::view::transform([](auto&& widths){return widths.first * widths.second; })
  //                     | ranges::view::transform([](auto&& square){return square / electron; });
  //                     | ranges::view::chunk(size);
  //
  return result;
}

