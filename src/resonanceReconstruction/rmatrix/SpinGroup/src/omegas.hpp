template < typename Range >
auto omegas( const Energy& energy, const Range& coulombShifts ) const {

  return ranges::view::zip_with(
             [] ( const double w, const double phi )
                { return std::exp( std::complex< double >( 0.0, w - phi ) ); },
             coulombShifts, this->phaseShifts( energy ) );
}
