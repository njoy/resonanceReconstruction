template< typename T >
constexpr auto root( T t ){
  using Root =
    decltype( pow( typename decltype( 1.0 * t )::Units(), Ratio<1,2> ) );
  Quantity< Root > r; r.value = 1.0;
  return r;
}
