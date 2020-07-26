private:

/**
 *  @brief Private constructor
 */
CompoundSystem( std::vector< unresolved::SpinGroup >&& groups,
                std::vector< Energy >&& energies,
                unsigned int lmax ) :
  groups_( std::move( groups ) ), energies_( std::move( energies ) ),
  lmax_( lmax ) {

    verifySpinGroups( this->groups_ );
  }


public:

/**
 *  @brief Constructor
 *
 *  @param[in] groups   all l,J groups that apply to the compound system
 */
CompoundSystem( std::vector< unresolved::SpinGroup >&& groups ) :
  CompoundSystem( std::move( groups ), generateGrid( groups ),
                  getLMax( groups ) ) {}
