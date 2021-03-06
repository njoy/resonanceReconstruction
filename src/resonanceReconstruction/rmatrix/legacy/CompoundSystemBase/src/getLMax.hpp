static
unsigned int getLMax( const std::vector< SpinGroupType >& groups ) {

  if ( groups.size() != 0 ) {

    return ranges::max(
             groups
               | ranges::view::transform(
                   [] ( const auto& spingroup )
                      { return spingroup.orbitalAngularMomentum(); } ) );
  }
  return 0;
}
