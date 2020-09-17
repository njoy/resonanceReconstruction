/**
 *  @class
 *  @brief The compound system
 *
 *  This class contains the unresolved resonance parameters for legacy ENDF
 *  data. It can be used to generate cross section values at any energy but it
 *  also provides the energy grid that should be used (this grid is generated
 *  on request by the user).
 */
class CompoundSystem : protected CompoundSystemBase< unresolved::SpinGroup > {

public:

  /* constructor */
  using CompoundSystemBase::CompoundSystemBase;

  /* methods */
  using CompoundSystemBase::spinGroups;
  using CompoundSystemBase::evaluate;

  #include "resonanceReconstruction/rmatrix/legacy/unresolved/CompoundSystem/src/grid.hpp"
};
