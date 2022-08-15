auto resonances() const {
  return
    Parent::resonances
    | ranges::views::transform( []( const auto& r ) -> const Resonance&
                                { return static_cast< const Resonance& >( r ); } );
}
