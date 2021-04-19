#ifndef NJOY_R2_RMATRIX_CHANNELRADII
#define NJOY_R2_RMATRIX_CHANNELRADII

// system includes
#include <variant>

// other includes
#include "utility/overload.hpp"
#include "resonanceReconstruction/Quantity.hpp"
#include "resonanceReconstruction/rmatrix/ChannelRadiusTable.hpp"

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {

/**
 *  @class
 *  @brief Channel radii used in wave function calculations
 *
 *  The penetrability P, shift factor S and phase shift phi require knowledge
 *  of the channel radius in their calculation. The ChannelRadii class provides
 *  these radii for each on of these.
 */
class ChannelRadii {

public:

  /* alias */
  using ChannelRadiusVariant = std::variant< ChannelRadius,
                                             ChannelRadiusTable >;

private:
  
  /* fields */
  ChannelRadiusVariant penetrability_;
  ChannelRadiusVariant shiftFactor_;
  ChannelRadiusVariant phaseShift_;

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/ChannelRadii/src/ctor.hpp"

  /**
   *  @brief Return the channel radius for the penetrability P
   *
   *  @param[in] energy   the energy for which the radius must be given
   */
  ChannelRadius penetrabilityRadius( const Energy& energy ) const {

    return std::visit(
             njoy::utility::overload{ [&] ( const ChannelRadius& radius )
                                          { return radius; },
                                      [&] ( const ChannelRadiusTable& table )
                                          { return table( energy ); } },
             this->penetrability_ );
  }

  /**
   *  @brief Return the channel radius for the shift factor S
   *
   *  @param[in] energy   the energy for which the radius must be given
   */
  ChannelRadius shiftFactorRadius( const Energy& energy ) const {

    return std::visit(
             njoy::utility::overload{ [&] ( const ChannelRadius& radius )
                                          { return radius; },
                                      [&] ( const ChannelRadiusTable& table )
                                          { return table( energy ); } },
             this->shiftFactor_ );
  }

  /**
   *  @brief Return the channel radius for the phase shift phi
   *
   *  @param[in] energy   the energy for which the radius must be given
   */
  ChannelRadius phaseShiftRadius( const Energy& energy ) const {

    return std::visit(
             njoy::utility::overload{ [&] ( const ChannelRadius& radius )
                                          { return radius; },
                                      [&] ( const ChannelRadiusTable& table )
                                          { return table( energy ); } },
             this->phaseShift_ );
  }
};

} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
