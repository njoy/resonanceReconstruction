#ifndef NJOY_R2_RMATRIX_RESONANCETABLE
#define NJOY_R2_RMATRIX_RESONANCETABLE

// system includes
#include <vector>

// other includes
#include "Log.hpp"
#include "range/v3/algorithm/for_each.hpp"
#include "range/v3/view/all.hpp"
#include "resonanceReconstruction/Quantity.hpp"
#include "resonanceReconstruction/rmatrix/Resonance.hpp"
#include "resonanceReconstruction/rmatrix/ChannelID.hpp"

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {

/**
 *  @class
 *  @brief Resonance parameters for a specific J,pi value
 */
class ResonanceTable {

  /* fields */
  std::vector< ChannelID > channels_;
  std::vector< Resonance > widths_;

  /* auxiliary functions */
  #include "resonanceReconstruction/rmatrix/ResonanceTable/src/verifyTable.hpp"

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/ResonanceTable/src/ctor.hpp"

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
  auto channels() const {

    return ranges::cpp20::views::all( this->channels_ );
  }

  /**
   *  @brief Return the resonances
   */
  auto resonances() const {

    return ranges::cpp20::views::all( this->widths_ );
  }

  /**
   *  @brief Return the resonance energies
   */
  auto energies() const {

    return this->resonances()
             | ranges::cpp20::views::transform(
                   [] ( const auto& resonance ) -> decltype(auto)
                      { return resonance.energy(); } );
  }
};

} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
