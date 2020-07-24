/**
 *  @brief Constructor
 *
 *  @param[in] incident   the incident channel data for this l,J pair
 *  @param[in] table      the table of resonance parameters for this l,J pair
 */
SpinGroup( Channel< Neutron >&& incident, unresolved::ResonanceTable&& table ) :
  incident_( std::move( incident ) ),
  parameters_( std::move( table ) ) {}
