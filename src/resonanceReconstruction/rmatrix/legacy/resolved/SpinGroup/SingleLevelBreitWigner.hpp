/**
 *  @class
 *  @brief Resolved resonance data for a specific l,J spin group using SLBW
 *
 *  This class contains the resolved resonance parameters and the associated
 *  incident channel for legacy SLBW ENDF data.
 */
template <>
class SpinGroup< SingleLevelBreitWigner > : protected SpinGroupBase< resolved::ResonanceTable > {

  /* fields */
  Energy qx_;

  std::array< ReactionID, 3 > reactions_;

protected:

  const ReactionID& elastic() const { return this->reactions_[0]; }
  const ReactionID& capture() const { return this->reactions_[1]; }
  const ReactionID& fission() const { return this->reactions_[2]; }

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/legacy/resolved/SpinGroup/SingleLevelBreitWigner/src/ctor.hpp"

  /* methods */

  using SpinGroupBase::incidentChannel;
  using SpinGroupBase::orbitalAngularMomentum;
  using SpinGroupBase::totalAngularMomentum;
  using SpinGroupBase::resonanceTable;

  /**
   *  @brief Return competitive Q value
   */
  const Energy& QX() const { return this->qx_; }

  #include "resonanceReconstruction/rmatrix/legacy/resolved/SpinGroup/SingleLevelBreitWigner/src/grid.hpp"
  #include "resonanceReconstruction/rmatrix/legacy/resolved/SpinGroup/SingleLevelBreitWigner/src/evaluate.hpp"
};
