/**
 *  @brief Constructor
 *
 *  @param[in] groups   the different spin groups in the compound system
 */
CompoundSystem( std::vector< SpinGroup< Formalism, BoundaryOption > >&& groups ) :
  groups_( std::move( groups ) ) {

  verifySpinGroups( this->groups_ );
}

