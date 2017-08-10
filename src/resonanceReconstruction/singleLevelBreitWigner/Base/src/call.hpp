template< typename T >
static bool isTemperature( T&& t ){
  return hana::is_valid( []( auto&& t )->
                         decltype( Quantity<Kelvin>{} + t ){} )(t);
}

template< typename Functor, typename... Functors >
auto operator()( Quantity<ElectronVolts> energy,
                 // from child call
                 Functor&& functor,
                 Functors&&... functors ) const ->
  std::enable_if_t< not isTemperature(functor), CrossSection >{
  const auto crossSections =
    this->lvalues
    | ranges::view::transform
      ( [&, psichi = this->psichi( energy ) ]( auto&& lvalue ){
        return lvalue( energy,
                       psichi,
                       std::forward< Functor >( functor ),
                       std::forward< Functors >( functors )... ); } );

  return ranges::accumulate( crossSections | ranges::view::drop_exactly(1),
                             crossSections.front() );
}

/*
template< typename... Functors >
CrossSection operator()( Quantity<ElectronVolts> energy,
                         Quantity<Kelvin> temperature,
                         // from child call
                         Functors&&... functors ) const {
    const auto crossSections =
    this->lvalues
    | ranges::view::transform
      ( [&, psichi = this->psichi( energy, temperature ) ]( auto&& lvalue ){
        return lvalue( energy,
                       psichi,
                       std::forward< Functor >( functor ),
                       std::forward< Functors >( functors )... ); } );

  return ranges::accumulate( crossSections | ranges::view::drop_exactly(1),
                             crossSections.front() );
}
*/
