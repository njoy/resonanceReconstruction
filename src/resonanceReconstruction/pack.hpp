template< typename... Ts >
struct Pack;

template< typename... Ts >
auto pack( Ts&&... ts ) -> Pack< std::decay_t<Ts>... >;

template< typename... Ts >
struct Pack {
  std::tuple<Ts...> data;

  Pack( Ts... data ) : data( std::make_tuple(data...) ){}

  static constexpr auto indices =
    std::make_index_sequence< std::tuple_size< decltype(data) >::value >{};

  template< typename... Us, std::size_t... indices >
  auto add( const Pack< Us... >& other, const std::index_sequence< indices... > ) const {
    static_assert( std::tuple_size< decltype( data ) >::value
                   == std::tuple_size< decltype( other.data ) >::value, "" );
    return pack( ( std::get<indices>(data)
                   + std::get<indices>(other.data) )...  );
  }

  template< typename... Us, std::size_t... indices >
  auto subtract( const Pack< Us... >& other, const std::index_sequence< indices... > ) const {
    static_assert( std::tuple_size< decltype( data ) >::value
                   == std::tuple_size< decltype( other.data ) >::value, "" );
    return pack( ( std::get<indices>(data)
                              - std::get<indices>(other.data) )...  );
  }

  template< typename Scalar, std::size_t... indices >
  auto multiply( const Scalar& scalar, const std::index_sequence< indices... > ) const {
    return pack( ( std::get<indices>(data) * scalar )...  );
  }

  template< typename... Us >
  auto operator+( const Pack< Us... >& other ) const {
    return add( other, indices );
  }

  template< typename... Us >
  auto operator-( const Pack< Us... >& other ) const {
    return subtract( other, indices );
  }

  template< typename Scalar >
  auto operator*( const Scalar& scalar ) const {
    return multiply( scalar, indices );
  }
};

template< typename... Ts >
auto pack( Ts&&... ts ) -> Pack< std::decay_t<Ts>... > {
  return Pack< std::decay_t<Ts>... >{ std::forward<Ts>(ts)... };
}
