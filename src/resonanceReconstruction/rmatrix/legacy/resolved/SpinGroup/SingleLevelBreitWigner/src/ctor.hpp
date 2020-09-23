/**
 *  @brief Constructor
 *
 *  @param[in] incident   the incident channel data for this l,J pair
 *  @param[in] table      the table of resonance parameters for this l,J pair
 */
SpinGroup( Channel< Neutron >&& incident, resolved::ResonanceTable&& table,
           const Energy& qx ) :
  SpinGroupBase( std::move( incident ), std::move( table ) ),
  qx_( qx ),
  elastic_( table.resonances().size() ),
  capture_( table.resonances().size() ),
  fission_( table.resonances().size() ),
  total_( table.resonances().size() ),
  delta_( table.resonances().size() ) {}
