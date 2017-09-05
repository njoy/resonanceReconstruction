template< typename T >
static auto isTemperatureImpl( T&& ){
  return hana::is_valid( []( auto&& t )->
                         decltype( Quantity<Kelvin>{} + t ){} );
}

template< typename T >
static constexpr bool isTemperature( T&& t ){
  return decltype( isTemperatureImpl(t) ){};
}
