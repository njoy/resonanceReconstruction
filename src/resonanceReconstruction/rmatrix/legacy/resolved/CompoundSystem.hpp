/**
 *  @class
 *  @brief The compound system for legacy resolved SLBW and MLBW
 *
 *  This class contains the resolved resonance parameters for legacy ENDF
 *  data using SLBW and MLBW. It can be used to generate cross section values
 *  at any energy but it also provides the minimal energy grid that should be
 *  used (this grid is generated on request by the user).
 */
template < typename Formalism >
class CompoundSystem : protected CompoundSystemBase< SpinGroup< Formalism > > {

public:

  /* constructor */
  using CompoundSystemBase< SpinGroup< Formalism > >::CompoundSystemBase;

  /* methods */
  using CompoundSystemBase< SpinGroup< Formalism > >::spinGroups;
  using CompoundSystemBase< SpinGroup< Formalism > >::evaluate;

  #include "resonanceReconstruction/rmatrix/legacy/resolved/CompoundSystem/src/grid.hpp"
};
