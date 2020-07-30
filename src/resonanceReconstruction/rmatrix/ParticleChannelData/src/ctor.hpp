/**
 *  @brief Constructor
 *
 *  @param[in] channel    the particle channel
 *  @param[in] energies   the resonanc energies
 *  @param[in] widths     the reduced resonance widths
 */
ParticleChannelData( const ParticleChannel& channel,
                     std::vector< Energy >&& energies,
                     std::vector< ReducedWidth >&& widths ) :
  channel_( channel ), energies_( std::move( energies ) ),
  widths_( std::move( widths ) ) {

    verifySize( this->energies_, this->widths_ );
  }
