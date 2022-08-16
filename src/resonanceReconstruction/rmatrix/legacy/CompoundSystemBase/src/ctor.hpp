private:

/**
 *  @brief Private constructor
 */
CompoundSystemBase( std::vector< SpinGroupType >&& groups,
                    std::vector< ReactionID >&& reactions,
                    unsigned int lmax ) :
  groups_( std::move( groups ) ), reactions_( std::move( reactions ) ),
  lmax_( lmax ) {

    verifySpinGroups( this->groups_ );
  }


public:

/**
 *  @brief Constructor
 *
 *  @param[in] groups   all l,J groups that apply to the compound system
 */
CompoundSystemBase( std::vector< SpinGroupType >&& groups ) :
  CompoundSystemBase( std::move( groups ),
                      makeReactionIDs( groups ),
                      getLMax( groups ) ) {}
