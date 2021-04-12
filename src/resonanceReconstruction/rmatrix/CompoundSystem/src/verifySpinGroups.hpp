static
void verifySpinGroups(
    const std::vector< SpinGroup< Formalism, BoundaryOption > >& groups ) {

  // verify that there is at least one group

  if ( groups.size() == 0 ) {

    Log::error( "The number of spin groups in a compound system cannot be 0." );
    throw std::exception();
  }

  // verify that each group is unique (i.e. each Jpi occurs only once)

  auto toHalfIntegerString =
      [] ( const double a )
         { double half;
           return std::modf( a, &half ) == 0. ?
                      std::to_string( static_cast<int>( half ) ) :
                      std::to_string( 2 * static_cast<int>( half ) + 1 )
                        + "/2"; };

  const auto getQuantumNumbers = [] ( const auto& group ) {
    return std::visit( [] ( const auto& entry )
                          { return entry.quantumNumbers(); },
                       group.channels().front() );
  };

  const auto numbers = groups | ranges::views::transform( getQuantumNumbers );

  const auto verifyUniqueSpinGroup = [&] ( const auto& reference ) {

    if ( ranges::cpp20::count_if(
             numbers,
             [&] ( const auto& current )
                 { return ( current.totalAngularMomentum() ==
                                reference.totalAngularMomentum() ) &&
                          ( current.parity() ==
                                reference.parity() ); } ) > 1 ) {

      Log::error( "Spin groups in the compound system do not seem to be "
                  "unique." );
      Log::info( "There are at least two spin groups with Jpi={}",
                 toHalfIntegerString( reference.totalAngularMomentum() )
                     + ( reference.parity() > 0 ? "+" : "-" ) );
      throw std::exception();
    }
  };

  ranges::for_each( numbers, verifyUniqueSpinGroup );
}
