template < typename Formalism > class SpinGroup;

/**
 *  @class
 *  @brief Unresolved resonance data for a specific l,J spin group
 *
 *  This class contains the unresolved resonance parameters and the associated
 *  incident channel for legacy ENDF data.
 */
template <>
class SpinGroup< SingleLevelBreitWigner > : protected SpinGroupBase< resolved::ResonanceTable > {

  /* fields */
  Energy qx_;

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/legacy/resolved/SpinGroup/src/ctor.hpp"

  /* methods */

  using SpinGroupBase::incidentChannel;
  using SpinGroupBase::orbitalAngularMomentum;
  using SpinGroupBase::totalAngularMomentum;
  using SpinGroupBase::resonanceTable;

  /**
   *  @brief Return competitive Q value
   */
  const Energy& QX() const { return this->qx_; }

  #include "resonanceReconstruction/rmatrix/legacy/resolved/SpinGroup/src/evaluate.hpp"
};
