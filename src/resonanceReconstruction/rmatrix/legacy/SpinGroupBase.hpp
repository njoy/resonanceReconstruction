/**
 *  @class
 *  @brief Resonance data for a specific l,J spin group
 */
template < typename ResonanceTableType > class SpinGroupBase {

  /* fields */
  Channel< Neutron > incident_;
  ResonanceTableType table_;

public:

  /* constructor */

  /**
   *  @brief Constructor
   *
   *  @param[in] incident   the incident channel data for this l,J pair
   *  @param[in] table      the table of resonance parameters for this l,J pair
   */
  SpinGroupBase( Channel< Neutron >&& incident, ResonanceTableType&& table ) :
    incident_( std::move( incident ) ),
    table_( std::move( table ) ) {}

  /* methods */

  /**
   *  @brief Return the incident channel data
   *
   *  REMARK: not all information in the channel will be correct (channel spin,
   *          target parity, target charge, etc. since this data is not
   *          required for the legacy resolved/unresolved resonances)
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
  const ResonanceTableType& resonanceTable() const { return this->table_; }
};
