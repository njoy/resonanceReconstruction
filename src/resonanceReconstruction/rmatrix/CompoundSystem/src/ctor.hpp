/**
 *  @brief Constructor
 *
 *  @param[in] groups   the different spin groups in the compound system
 */
CompoundSystem( std::vector< SpinGroup< Formalism, BoundaryOption > >&& groups ) :
  groups_( std::move( groups ) ) {

  verifySpinGroups( this->groups_ );
}

/**
 *  @brief Constructor
 *
 *  @param[in] channels   the different channels and their widths that make up
 *                        the compound system
 */
CompoundSystem( std::vector< ParticleChannelData >&& channels ) :
  CompoundSystem( makeSpinGroups( std::move( channels ) ) ) {}
