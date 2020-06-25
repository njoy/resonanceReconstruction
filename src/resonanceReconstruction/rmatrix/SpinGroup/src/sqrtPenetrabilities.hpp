template < typename Range >
auto sqrtPenetrabilities( const Range& penetrabilities ) const {

  return penetrabilities
             | ranges::view::transform(
                   [] ( const double penetrability ) -> double
                      { return std::sqrt( penetrability ); } );
}
