private:

/**
 *  @brief Private constructor
 */
CompoundSystem( std::vector< unresolved::SpinGroup >&& groups,
                unsigned int lmax ) :
  groups_( std::move( groups ) ), lmax_( lmax ) {}


public:

/**
 *  @brief Constructor
 *
 *  @param[in] groups   all l,J groups that apply to the compound system
 */
CompoundSystem( std::vector< unresolved::SpinGroup >&& groups ) :
  CompoundSystem( std::move( groups ), getLMax( groups ) ) {}
