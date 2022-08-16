/**
 *  @class
 *  @brief Resolved resonance data for a specific l,J spin group using MLBW
 *
 *  This class contains the resolved resonance parameters and the associated
 *  incident channel for legacy MLBW ENDF data.
 */
template <>
class SpinGroup< MultiLevelBreitWigner > : protected SpinGroup< SingleLevelBreitWigner > {

public:

  /* constructor */
  using SpinGroup< SingleLevelBreitWigner >::SpinGroup;

  /* methods */
  using SpinGroup< SingleLevelBreitWigner >::incidentChannel;
  using SpinGroup< SingleLevelBreitWigner >::orbitalAngularMomentum;
  using SpinGroup< SingleLevelBreitWigner >::totalAngularMomentum;
  using SpinGroup< SingleLevelBreitWigner >::resonanceTable;
  using SpinGroup< SingleLevelBreitWigner >::QX;
  using SpinGroup< SingleLevelBreitWigner >::grid;
  using SpinGroup< SingleLevelBreitWigner >::reactionIDs;
  using SpinGroup< SingleLevelBreitWigner >::hasFission;

  #include "resonanceReconstruction/rmatrix/legacy/resolved/SpinGroup/MultiLevelBreitWigner/src/evaluate.hpp"
};
