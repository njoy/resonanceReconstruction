/**
 *  @class
 *  @brief Resonance parameters for a specific J,pi value
 *
 *  The ResonanceTable is used to store resonance data and extract specific
 *  data for use in resonance reconstruction
 */
class ResonanceTable {

  /* fields */
  std::vector< ChannelID > channels_;
  std::vector< Energy > energies_;
  std::vector< std::vector< ReducedWidth > > widths_;

  //! @todo store index for a channel at construction of the table

public:

  /* constructor */

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
                  std::vector< Energy >&& energies,
                  std::vector< std::vector< ReducedWidth > >&& widths ) :
      channels_( std::move( channels ) ),
      energies_( std::move( energies ) ),
      widths_( std::move( widths ) ) {}

  /**
   *  @brief Return the resonance energies
   */
  auto energies() const {
    return ranges::make_iterator_range( this->energies_.begin(),
                                        this->energies_.end() );
  }

  /**
   *  @brief Return the reduced widths for a specific channel
   *
   *  @param[in] channel   the channel ID
   */
  auto reducedWidths( const ChannelID& channel ) const {
    unsigned int index = std::distance( this->channels_.begin(),
                                        std::find( this->channels_.begin(),
                                                   this->channels_.end(),
                                                   channel ) );
    if ( index != this->channels_.size() ) {
      return this->widths_ 
               | ranges::view::transform( [&index] ( const auto& array )
                                                   { return array[index]; } );
    }
    else {
      throw std::exception();
    }
  }
};
