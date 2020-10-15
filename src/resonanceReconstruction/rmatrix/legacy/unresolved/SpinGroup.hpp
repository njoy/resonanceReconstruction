/**
 *  @class
 *  @brief Unresolved resonance data for a specific l,J spin group
 *
 *  This class contains the unresolved resonance parameters and the associated
 *  incident channel for legacy ENDF data.
 */
class SpinGroup : protected SpinGroupBase< unresolved::ResonanceTable > {

public:

  /* constructor */

  using SpinGroupBase::SpinGroupBase;

  /* methods */

  using SpinGroupBase::incidentChannel;
  using SpinGroupBase::orbitalAngularMomentum;
  using SpinGroupBase::totalAngularMomentum;
  using SpinGroupBase::resonanceTable;

  #include "resonanceReconstruction/rmatrix/legacy/unresolved/SpinGroup/src/evaluate.hpp"
};
