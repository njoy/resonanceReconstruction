/**
 *  @class
 *  @brief The base interface for a legacy compound system
 *
 *  This class contains the resolved or unresolved resonance parameters
 *  for legacy ENDF data.
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
