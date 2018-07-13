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
  std::vector< Resonance > widths_;

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
                  std::vector< Resonance >&& widths ) :
      channels_( std::move( channels ) ),
      widths_( std::move( widths ) ) {}

  /**
   *  @brief Return the number of channels
   */
  unsigned int numberChannels() const {
    return this->channels_.size();
  }

  /**
   *  @brief Return the number of resonances
   */
  unsigned int numberResonances() const {
    return this->widths_.size();
  }

  /**
   *  @brief Return the channel IDs
   */
  auto channels() const {
    return ranges::make_iterator_range( this->channels_.begin(),
                                        this->channels_.end() );
  }

  /**
   *  @brief Return the resonances
   */
  auto resonances() const {
    return ranges::view::all( this->widths_ );
  }

  /**
   *  @brief Return the resonance energies
   */
  auto energies() const {
    return this->resonances()
             | ranges::view::transform( [] ( const auto& resonance )
                                           { return resonance.energy(); } );
  }
};
