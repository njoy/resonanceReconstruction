/**
 *  @brief Constructor
 *
 *  The ResonanceTable class takes the reduced widths for all channels in a
 *  spin group. As the channels can differe from case to case, channels are
 *  identified using their channel IDs. Data can be extracted from the table
 *  using these IDs.
 *
 *  @param[in] channels   the channel IDs (nc values)
 *  @param[in] energies   the energies for which resonances are defined (ne
 *                        values)
 *  @param[in] widths     the reduced widths for each energy and channel
 *                        (ne arrays of nc values)
 */
ResonanceTable( std::vector< ChannelID >&& channels,
                std::vector< Resonance >&& widths ) :
    channels_( std::move( channels ) ),
    widths_( std::move( widths ) ) {

    verifyTable( this->channels_, this->widths_ );
}

