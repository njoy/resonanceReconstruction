/**
 *  @brief Constructor
 *
 *  @param[in] incident   the incident channel data for this l,J pair
 *  @param[in] table      the table of resonance parameters for this l,J pair
 */
SpinGroup( Channel< Neutron >&& incident, ResonanceTableType&& table,
           const Energy& qx ) :
  SpinGroupBase( std::move( incident ), std::move( table ) ),
  qx_( qx ) {}
