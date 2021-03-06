static
void verifySpinGroups( const std::vector< SpinGroupType >& groups ) {

  // verify that there is at least one group

  if ( groups.size() == 0 ) {

    Log::error( "The number of spin groups in a compound system cannot be 0." );
    throw std::exception();
  }

  // verify that each group is unique (i.e. each l,J occurs only once)

  auto toHalfIntegerString =
      [] ( const double a )
         { double half;
           return std::modf( a, &half ) == 0. ?
                      std::to_string( static_cast<int>( half ) ) :
                      std::to_string( 2 * static_cast<int>( half ) + 1 )
                        + "/2"; };

  const auto getQuantumNumbers = [] ( const auto& group ) {

    return group.incidentChannel().quantumNumbers();
  };

  const auto numbers = groups | ranges::view::transform( getQuantumNumbers );

  const auto verifyUniqueSpinGroup = [&] ( const auto& reference ) {

    if ( ranges::count_if(
             numbers,
             [&] ( const auto& current )
                 { return ( current.orbitalAngularMomentum() ==
                                reference.orbitalAngularMomentum() ) &&
                          ( current.totalAngularMomentum() ==
                                reference.totalAngularMomentum() ); } ) > 1 ) {

      Log::error( "The spin groups in the compound system do not "
                  "seem to be unique." );
      Log::info( "There are at least two spin groups with l,J={},{}",
                 reference.orbitalAngularMomentum(),
                 toHalfIntegerString( reference.totalAngularMomentum() ) );
      throw std::exception();
    }
  };

  ranges::for_each( numbers, verifyUniqueSpinGroup );
}
