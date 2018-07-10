/**
 *  @class
 *  @brief Channel radii used in wave function calculations
 *
 *  The penetrability P, shift factor S and phase shift phi require knowledge
 *  of the channel radius in their calculation. The ChannelRadii class provides
 *  these radii for each on of these.
 */
class ChannelRadii {

  //! @todo make these energy dependent

  /* fields */
  ChannelRadius penetrability_;
  ChannelRadius shiftFactor_;
  ChannelRadius phaseShift_;

public:

  /* constructor */
  ChannelRadii( const ChannelRadius& radius ) :
      penetrability_( radius ),
      shiftFactor_( radius ),
      phaseShift_( radius ) {}

  ChannelRadii( const ChannelRadius& trueRadius,
                const ChannelRadius& effectiveRadius ) :
      penetrability_( trueRadius ),
      shiftFactor_( trueRadius ),
      phaseShift_( effectiveRadius ) {}

  ChannelRadii( const ChannelRadius& penetrationRadius,
                const ChannelRadius& shiftFactorRadius,
                const ChannelRadius& phaseShiftRadius ) :
      penetrability_( penetrationRadius ),
      shiftFactor_( shiftFactorRadius ),
      phaseShift_( phaseShiftRadius ) {}

  /**
   *  @brief Return the channel radius for the penetrability P
   */
  const ChannelRadius& penetrabilityRadius( /*const Energy& energy*/ ) const {
    return this->penetrability_;
  }

  /**
   *  @brief Return the channel radius for the shift factor S
   */
  const ChannelRadius& shiftFactorRadius( /*const Energy& energy*/ ) const {
    return this->shiftFactor_;
  }

  /**
   *  @brief Return the channel radius for the phase shift phi
   */
  const ChannelRadius& phaseShiftRadius( /*const Energy& energy*/ ) const {
    return this->phaseShift_;
  }

};
