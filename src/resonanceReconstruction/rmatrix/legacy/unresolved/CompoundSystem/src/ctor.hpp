/**
 *  @brief Constructor
 *
 *  @param[in] groups          all l,J groups that apply to the compound system
 *  @param[in] interpolation   interpolation scheme to be used for cross
 *                             sections
 */
CompoundSystem( std::vector< unresolved::SpinGroup >&& groups,
                int interpolation ) :
  CompoundSystemBase( std::move( groups ) ), interpolation_( interpolation ) {}
