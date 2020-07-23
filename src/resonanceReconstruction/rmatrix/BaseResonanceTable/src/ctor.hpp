/**
 *  @brief Constructor
 *
 *  @param[in] channels   the channel IDs (nc values)
 *  @param[in] energies   the energies for which resonances are defined (ne
 *                        values)
 *  @param[in] widths     the reduced widths for each energy and channel
 *                        (ne arrays of nc values)
 */
BaseResonanceTable( std::vector< ChannelID >&& channels,
                    std::vector< ResonanceType >&& widths ) :
    channels_( std::move( channels ) ),
    widths_( std::move( widths ) ) {

    verifyTable( this->channels_, this->widths_ );
}
