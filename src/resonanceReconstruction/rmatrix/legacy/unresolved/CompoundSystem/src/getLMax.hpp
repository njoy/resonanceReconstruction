static
unsigned int getLMax( const std::vector< SpinGroup >& groups ) {

  return ranges::max(
           groups
             | ranges::view::transform(
                 [] ( const auto& spingroup )
                    { return spingroup.orbitalAngularMomentum(); } ) );
}
