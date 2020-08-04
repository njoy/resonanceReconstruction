/**
 *  @class
 *  @brief The compound nucleus system
 */
template < typename Formalism, typename BoundaryOption >
class CompoundSystem {

  /* fields */
  std::vector< SpinGroup< Formalism, BoundaryOption > > groups_;

  /* auxiliary functions */
  #include "resonanceReconstruction/rmatrix/CompoundSystem/src/makeSpinGroups.hpp"
  #include "resonanceReconstruction/rmatrix/CompoundSystem/src/verifySpinGroups.hpp"

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/CompoundSystem/src/ctor.hpp"

  auto spinGroups() const { return ranges::view::all( this->groups_ ); }

  //#include "resonanceReconstruction/rmatrix/CompoundSystem/src/switchIncidentPair.hpp"
  #include "resonanceReconstruction/rmatrix/CompoundSystem/src/evaluate.hpp"
  #include "resonanceReconstruction/rmatrix/CompoundSystem/src/evaluateTMatrix.hpp"
  #include "resonanceReconstruction/rmatrix/CompoundSystem/src/grid.hpp"
};
