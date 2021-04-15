template < typename Range >
auto omegas( const Energy& energy, const Range& coulombShifts ) const {

  const auto omega = [] ( const double w, const double phi ) {

    return std::exp( std::complex< double >( 0.0, w - phi ) );
  };

  return ranges::views::zip_with( omega,
                                  coulombShifts, this->phaseShifts( energy ) );
}
