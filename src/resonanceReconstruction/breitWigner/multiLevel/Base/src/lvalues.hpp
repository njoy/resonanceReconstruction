auto lvalues() const {
  return
    Parent::lvalues
    | ranges::view::transform( []( const auto& l ) -> const Lvalue&
                               { return static_cast< const Lvalue& >( l ); } );
}
