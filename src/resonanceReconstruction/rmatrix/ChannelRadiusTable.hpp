#ifndef NJOY_R2_RMATRIX_CHANNELRADIUSTABLE
#define NJOY_R2_RMATRIX_CHANNELRADIUSTABLE

// system includes
#include <memory>

// other includes
#include "resonanceReconstruction/Quantity.hpp"
#include "resonanceReconstruction/rmatrix/Table.hpp"

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {

/**
 *  @class
 *  @brief Energy dependent channel radius table
 *
 *  This class contains a shared pointer to a MultiRegionTable because the
 *  energy dependent scattering radius can be shared by multiple channels. ENDF
 *  only defines one energy dependent channel radius which is then used in all
 *  lJ spin groups for SLBW, MLBW, RM resolved resonances and unresolved
 *  resonances. GNDS is more broader than this and allows for energy dependent
 *  channel radii to be defined up to the spin group level.
 */
class ChannelRadiusTable {

  /* aliases */
  using RadiusTable = MultiRegionTable< Energy, ChannelRadius >;

  /* fields */
  std::shared_ptr< RadiusTable > table_;

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/ChannelRadiusTable/src/ctor.hpp"

  /**
   *  @brief Return the channel radius at the requested energy
   *
   *  @param[in] energy   the energy for which the radius must be given
   */
  ChannelRadius operator()( const Energy& energy ) const {

    return this->table_->operator()( energy );
  }
};

} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
