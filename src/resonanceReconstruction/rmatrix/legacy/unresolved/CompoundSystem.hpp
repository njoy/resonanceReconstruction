/**
 *  @class
 *  @brief The compound system
 *
 *  This class contains the unresolved resonance parameters for legacy ENDF
 *  data.
 */
class CompoundSystem {

  /* fields */
  std::vector< SpinGroup > groups_;
  unsigned int lmax_;

  /* auxiliary functions */
  #include "resonanceReconstruction/rmatrix/legacy/unresolved/CompoundSystem/src/getLMax.hpp"

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/legacy/unresolved/CompoundSystem/src/ctor.hpp"

  /**
   *  @brief Return the l,J data
   */
  auto spinGroups() const {

    return ranges::view::all( this->groups_ );
  }

  #include "resonanceReconstruction/rmatrix/legacy/unresolved/CompoundSystem/src/evaluate.hpp"
};
