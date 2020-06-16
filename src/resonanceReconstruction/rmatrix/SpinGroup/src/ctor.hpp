private:

/**
 *  @brief Constructor
 *
 *  @param[in] incidentChannels   the indices of the incident channels
 *  @param[in] channels           the channels involved in the spin group
 *  @param[in] table              the table of resonance parameters
 */
SpinGroup( std::vector< unsigned int >&& incidentChannels,
           std::vector< ParticleChannel >&& channels,
           ResonanceTable&& table ) :
  rlmatrix_( table ),
  reactions_( makeReactionIdentifiers( channels,
                                       incidentChannels,
                                       Formalism() ) ),
  incident_( std::move( incidentChannels ) ),
  channels_( std::move( channels ) ),
  parameters_( std::move( table ) ) {

  verifyChannels( this->channels_ );
  verifyIncidentChannels( this->incident_ );
  verifyResonanceChannels( this->channels_, this->parameters_ );
}

public:

/**
 *  @brief Constructor
 *
 *  @param[in] incident   the incident particle pair
 *  @param[in] channels   the channels involved in the spin group
 *  @param[in] table      the table of resonance parameters
 */
SpinGroup( const ParticlePair& incident,
           std::vector< ParticleChannel >&& channels,
           ResonanceTable&& table ) :
  SpinGroup( determineIncidentChannels( incident.pairID(), channels ),
             std::move( channels ),
             std::move( table ) ) {}
