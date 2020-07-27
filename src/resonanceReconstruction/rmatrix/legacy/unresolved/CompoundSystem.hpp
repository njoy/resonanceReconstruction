/**
 *  @class
 *  @brief The compound system
 *
 *  This class contains the unresolved resonance parameters for legacy ENDF
 *  data. It can be used to generate cross section values at any energy but it
 *  also provides the energy grid that should be used (this grid is generated
 *  at construction time and can be retrieved by the user).
 */
class CompoundSystem {

  /* fields */
  std::vector< SpinGroup > groups_;
  std::vector< Energy > energies_;
  unsigned int lmax_;

  /* auxiliary functions */
  #include "resonanceReconstruction/rmatrix/legacy/unresolved/CompoundSystem/src/getLMax.hpp"
  #include "resonanceReconstruction/rmatrix/legacy/unresolved/CompoundSystem/src/verifySpinGroups.hpp"

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/legacy/unresolved/CompoundSystem/src/ctor.hpp"

  /**
   *  @brief Return the l,J data
   */
  auto spinGroups() const { return ranges::view::all( this->groups_ ); }

  #include "resonanceReconstruction/rmatrix/legacy/unresolved/CompoundSystem/src/evaluate.hpp"
  #include "resonanceReconstruction/rmatrix/legacy/unresolved/CompoundSystem/src/grid.hpp"
};
