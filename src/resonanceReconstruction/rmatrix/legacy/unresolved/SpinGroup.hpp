/**
 *  @class
 *  @brief Unresolved resonance data for a specific l,J spin group
 *
 *  This class contains the unresolved resonance parameters and the associated
 *  incident channel for legacy ENDF data.
 */
class SpinGroup {

  /* fields */
  Channel< Neutron > incident_;
  unresolved::ResonanceTable parameters_;

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/legacy/unresolved/SpinGroup/src/ctor.hpp"

  /**
   *  @brief Return the incident channel data
   */
  const Channel< Neutron >& incidentChannel() const { return this->incident_; }

  /**
   *  @brief Return the orbital angular momentum l for this l,J pair
   */
  const OrbitalAngularMomentum& orbitalAngularMomentum() const {

    return this->incidentChannel().quantumNumbers().orbitalAngularMomentum();
  }

  /**
   *  @brief Return the total angular momentum J for this l,J pair
   */
  const TotalAngularMomentum& totalAngularMomentum() const {

    return this->incidentChannel().quantumNumbers().totalAngularMomentum();
  }

  /**
   *  @brief Return the resonance table
   */
  const unresolved::ResonanceTable& resonanceTable() const {

    return this->parameters_;
  }

  #include "resonanceReconstruction/rmatrix/legacy/unresolved/SpinGroup/src/evaluate.hpp"
};
