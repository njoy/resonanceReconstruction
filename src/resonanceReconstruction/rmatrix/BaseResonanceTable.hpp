/**
 *  @class
 *  @brief Resonance parameters for a specific J,pi or l,J value
 */
template < typename ResonanceType >
class BaseResonanceTable {

  /* fields */
  std::vector< ChannelID > channels_;
  std::vector< ResonanceType > widths_;

  /* auxiliary functions */
  #include "resonanceReconstruction/rmatrix/BaseResonanceTable/src/verifyTable.hpp"

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/BaseResonanceTable/src/ctor.hpp"

  /**
   *  @brief Return the number of channels
   */
  unsigned int numberChannels() const { return this->channels_.size(); }

  /**
   *  @brief Return the number of resonances
   */
  unsigned int numberResonances() const { return this->widths_.size(); }

  /**
   *  @brief Return the channel IDs
   */
  auto channels() const { return ranges::view::all( this->channels_ ); }

  /**
   *  @brief Return the resonances
   */
  auto resonances() const { return ranges::view::all( this->widths_ ); }

  /**
   *  @brief Return the resonance energies
   */
  auto energies() const {

    return this->resonances()
             | ranges::view::transform( [] ( const auto& resonance )
                                           { return resonance.energy(); } );
  }
};
