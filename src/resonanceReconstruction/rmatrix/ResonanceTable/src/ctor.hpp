/**
 *  @brief Constructor
 *
 *  The ResonanceTable class takes the reduced widths for all channels in a
 *  spin group. As the channels can differe from case to case, channels are
 *  identified using their channel IDs. Data can be extracted from the table
 *  using these IDs.
 *
 *  @param[in] channels     the channel identifiers
 *  @param[in] resonances   the resolved resonance parameters
 */
ResonanceTable( std::vector< ChannelID >&& channels,
                std::vector< Resonance >&& resonances ) :
    channels_( std::move( channels ) ),
    widths_( std::move( resonances ) ) {

    verifyTable( this->channels_, this->widths_ );
}
