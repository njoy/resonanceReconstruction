static
Matrix< std::complex< double > >
makeTemporaryMatrix( const ResonanceTable& table,
                     ReichMoore ) {

  // the temporary matrix is the R-matrix for the Reich-Moore approximation
  return Matrix< std::complex< double > >( table.numberChannels(),
                                           table.numberChannels() );
}

static
Matrix< std::complex< double > >
makeTemporaryMatrix( const ResonanceTable& table,
                     GeneralRMatrix ) {

  // the temporary matrix is a rectangular matrix with the widths of each
  // channel and resonance level

  Matrix< std::complex< double > > matrix( table.numberChannels(),
                                           table.numberResonances() );

  auto resonances = table.resonances();
  for ( unsigned int lambda = 0; lambda < table.numberResonances(); ++lambda ) {

    auto widths = resonances[ lambda ].widths();
    for ( unsigned int c = 0; lambda < table.numberChannels(); ++c ) {

      matrix( c, lambda ) = widths[c];
    }
  }

  return matrix;
}
