/**
 *  @class
 *  @brief Resonance parameters for a specific J,pi value
 */
class ResonanceTable {

  /* fields */
  std::vector< ChannelID > channels_;
  std::vector< Energy > energies__;
  std::vector< std::vector< ResonanceWidth > > widths__;

public:

  /* constructor */

/*  auto channelWidths( const ChannelID& channel ) const {
    unsigned int stride = std::distance( this->channels_.begin(),
                                         std::find( this->channels_.begin(),
                                                    this->channels_.end(),
                                                    channel ) );
    if ( stride != this->channels_.size() ) {
        return ;
    }
    else {
      throw std::exception();
    }
  }*/
};
