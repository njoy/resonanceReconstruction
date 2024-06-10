static
Matrix< std::complex< double > >
makeGMatrix( const ResonanceTable& table ) {

  Matrix< std::complex< double > > matrix( table.numberChannels(),
                                           table.numberResonances() );

  auto resonances = table.resonances();
  for ( unsigned int lambda = 0; lambda < table.numberResonances(); ++lambda ) {

    auto widths = resonances[ lambda ].widths();
    for ( unsigned int c = 0; c < table.numberChannels(); ++c ) {

      matrix( c, lambda ) = widths[c] / rootElectronVolt;
    }
  }

  return matrix;
}
