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
  reactions_( makeReactionIdentifiers( channels,
                                       incidentChannels.front() ) ),
  uMatrix_( channels.size(), channels.size() ),
  diagonalLMatrix_( channels.size() ),
  incident_( std::move( incidentChannels ) ),
  channels_( std::move( channels ) ),
  parameters_( std::move( table ) ) {

  //! @todo check if there are incident channels (i.e. size != 0)
  //! @todo check all ChannelQuantumNumbers in each Channel for equal J,pi
  //! @todo check if the incident channels have the same particle pair
  //! @todo check if number of widths is the same as the number of channels
  //! @todo check if channel IDs in the table and the channels correspond
}

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

