static
void verifySpinGroups(
    const std::vector< SpinGroup< Formalism, BoundaryOption > >& groups ) {

  // verify that there is at least one group

  if ( groups.size() == 0 ) {

    Log::error( "The number of spin groups in a compound system cannot be 0." );
    throw std::exception();
  }

  // verify that each group is unique (i.e. each Jpi occurs only once)

  const auto toHalfIntegerString = [] ( const double a ) {

    double half;
    return std::modf( a, &half ) == 0. ?
               std::to_string( static_cast<int>( half ) ) :
               std::to_string( 2 * static_cast<int>( half ) + 1 ) + "/2";
  };

  const auto getNumbers = [] ( const auto& group ) -> decltype(auto) {

    return std::visit( [] ( const auto& entry ) -> decltype(auto)
                          { return entry.quantumNumbers(); },
                       group.channels().front() );
  };

  const auto numbers = groups | ranges::cpp20::views::transform( getNumbers );

  const auto verifyUniqueSpinGroup = [&] ( const auto& reference ) {

    const auto equalToReferenceNumbers = [&] ( const auto& current ) {

      return ( current.totalAngularMomentum()
                  == reference.totalAngularMomentum() ) &&
             ( current.parity() == reference.parity() );
    };

    if ( ranges::cpp20::count_if( numbers, equalToReferenceNumbers ) > 1 ) {

      Log::error( "Spin groups in the compound system do not seem to be "
                  "unique." );
      Log::info( "There are at least two spin groups with Jpi={}",
                 toHalfIntegerString( reference.totalAngularMomentum() )
                     + ( reference.parity() > 0 ? "+" : "-" ) );
      throw std::exception();
    }
  };

  ranges::cpp20::for_each( numbers, verifyUniqueSpinGroup );
}
