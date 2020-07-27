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

  //! @todo make these energy dependent using a variant

  /* fields */
  ChannelRadius penetrability_;
  ChannelRadius shiftFactor_;
  ChannelRadius phaseShift_;

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/ChannelRadii/src/ctor.hpp"

  /**
   *  @brief Return the channel radius for the penetrability P
   *
   *  @param[in] energy   the energy for which the radius must be given (unused)
   */
  const ChannelRadius& penetrabilityRadius( const Energy& ) const {

    return this->penetrability_;
  }

  /**
   *  @brief Return the channel radius for the shift factor S
   *
   *  @param[in] energy   the energy for which the radius must be given (unused)
   */
  const ChannelRadius& shiftFactorRadius( const Energy& ) const {

    return this->shiftFactor_;
  }

  /**
   *  @brief Return the channel radius for the phase shift phi
   *
   *  @param[in] energy   the energy for which the radius must be given (unused)
   */
  const ChannelRadius& phaseShiftRadius( const Energy& ) const {

    return this->phaseShift_;
  }
};
