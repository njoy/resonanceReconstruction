/**
 *  @class
 *  @brief The compound nucleus system
 */
template < typename BoundaryOption >
class CompoundNucleus {

  /* fields */
  std::vector< SpinGroup< BoundaryOption > > groups_;

public:

  /* constructor */
  CompoundNucleus( std::vector< SpinGroup< BoundaryOption > >&& groups ) :
    groups_( std::move( groups ) ) {}

  auto spinGroups() const { return ranges::view::all( this->groups_ ); }

  void evaluate( const Energy& energy,
                 tsl::hopscotch_map< ReactionID, Quantity< Barn > >& result ) {

    ranges::for_each( this->groups_,
                      [&] ( auto& group )
                          { group.evaluate( energy, result ); } );
  }
};
