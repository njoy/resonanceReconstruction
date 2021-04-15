template < typename Range >
auto sqrtPenetrabilities( const Range& penetrabilities ) const {

  const auto sqrt = [] ( const double penetrability ) -> double {

    return std::sqrt( penetrability );
  };

  return penetrabilities | ranges::cpp20::views::transform( sqrt );
}
