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
  reactions_(
    [] ( const auto& channel ) -> std::array< ReactionID, 3 > {

      auto incident = channel.particlePair().particle().particleID();
      auto target = channel.particlePair().residual().particleID();
      return {{ ReactionID{ incident, target, ReactionType( "elastic" ) },
                ReactionID{ incident, target, ReactionType( "capture" ) },
                ReactionID{ incident, target, ReactionType( "fission" ) } }};
    }( incident ) ),
  elastic_( table.resonances().size() ),
  capture_( table.resonances().size() ),
  fission_( table.resonances().size() ),
  total_( table.resonances().size() ),
  delta_( table.resonances().size() ) {}
