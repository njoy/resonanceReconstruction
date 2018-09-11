/**
 *  @class
 *  @brief The compound nucleus system
 */
template < typename Formalism, typename BoundaryOption >
class CompoundSystem {

  /* fields */
  std::vector< SpinGroup< Formalism, BoundaryOption > > groups_;

  /* auxiliary functions */
  #include "resonanceReconstruction/rmatrix/CompoundSystem/src/verifySpinGroups.hpp"

public:

  /* constructor */
  CompoundSystem( std::vector< SpinGroup< Formalism, BoundaryOption > >&& groups ) :
    groups_( std::move( groups ) ) {

    //! @todo check for potential duplicate J,pi?
  }

  auto spinGroups() const { return ranges::view::all( this->groups_ ); }

  #include "resonanceReconstruction/rmatrix/CompoundSystem/src/switchIncidentPair.hpp"
  #include "resonanceReconstruction/rmatrix/CompoundSystem/src/evaluate.hpp"
};
