/**
 *  @brief Constructor
 *
 *  @param[in] channels   the channels involved in the spin group
 *  @param[in] table      the table of resonance parameters
 */
SpinGroup( std::vector< ParticleChannel >&& channels,
           ResonanceTable&& table ) :
  rlmatrix_( table ),
  reactions_( makeReactionIdentifiers( channels,
                                       Formalism() ) ),
  incident_( determineIncidentChannels( channels ) ),
  channels_( std::move( channels ) ),
  parameters_( std::move( table ) ) {

  verifyChannels( this->channels_ );
  verifyIncidentChannels( this->incident_ );
  verifyResonanceChannels( this->channels_, this->parameters_ );
}
