/**
 *  @class
 *  @brief Channel radii used in wave function calculations
 *
 *  The penetrability P, shift factor S and phase shift phi require knowledge
 *  of the channel radius in their calculation. The ChannelRadii class provides
 *  these radii for each on of these.
 *
 *  This class only accepts energy independent radii, but the interface does
 *  already provides a semblance of energy dependence when retrieving the
 *  radii.
 */
class ChannelRadii {

  /* alias */
  using ChannelRadiusVariant = std::variant< ChannelRadius,
                                             ChannelRadiusTable >;

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
             overload{ [&] ( const ChannelRadius& radius )
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
             overload{ [&] ( const ChannelRadius& radius )
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
             overload{ [&] ( const ChannelRadius& radius )
                           { return radius; },
                       [&] ( const ChannelRadiusTable& table )
                           { return table( energy ); } },
             this->phaseShift_ );
  }
};
