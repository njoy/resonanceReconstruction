/**
 *  @class
 *  @brief The compound system
 *
 *  This class contains the unresolved resonance parameters for legacy ENDF
 *  data. It can be used to generate cross section values at any energy but it
 *  also provides the energy grid that should be used (this grid is generated
 *  on request by the user).
 */
template < typename SpinGroupType >
class CompoundSystemBase {

  /* fields */
  std::vector< SpinGroupType > groups_;
  unsigned int lmax_;

  /* auxiliary functions */
  #include "resonanceReconstruction/rmatrix/legacy/CompoundSystemBase/src/getLMax.hpp"
  #include "resonanceReconstruction/rmatrix/legacy/CompoundSystemBase/src/verifySpinGroups.hpp"

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/legacy/CompoundSystemBase/src/ctor.hpp"

  /**
   *  @brief Return the l,J data
   */
  auto spinGroups() const { return ranges::view::all( this->groups_ ); }

  #include "resonanceReconstruction/rmatrix/legacy/CompoundSystemBase/src/evaluate.hpp"
};
