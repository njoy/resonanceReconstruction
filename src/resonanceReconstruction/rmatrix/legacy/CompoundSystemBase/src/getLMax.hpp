static
unsigned int getLMax( const std::vector< SpinGroupType >& groups ) {

  if ( groups.size() != 0 ) {

    return ranges::cpp20::max(
             groups
               | ranges::views::transform(
                   [] ( const auto& spingroup )
                      { return spingroup.orbitalAngularMomentum(); } ) );
  }
  return 0;
}
