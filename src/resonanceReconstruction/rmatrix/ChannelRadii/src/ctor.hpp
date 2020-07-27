/**
 *  @brief Constructor
 *
 *  @param[in] radius   the channel radius to be used for P, S and phi
 */
ChannelRadii( const ChannelRadius& radius ) :
    penetrability_( radius ),
    shiftFactor_( radius ),
    phaseShift_( radius ) {}

/**
 *  @brief Constructor
 *
 *  @param[in] trueRadius        the channel radius to be used for P and S
 *  @param[in] effectiveRadius   the channel radius to be used for phi
 */
ChannelRadii( const ChannelRadius& trueRadius,
             const ChannelRadius& effectiveRadius ) :
    penetrability_( trueRadius ),
    shiftFactor_( trueRadius ),
    phaseShift_( effectiveRadius ) {}

/**
 *  @brief Constructor
 *
 *  @param[in] penetrationRadius   the channel radius to be used for P
 *  @param[in] shiftFactorRadius   the channel radius to be used for S
 *  @param[in] phaseShiftRadius    the channel radius to be used for phi
 */
ChannelRadii( const ChannelRadius& penetrationRadius,
             const ChannelRadius& shiftFactorRadius,
             const ChannelRadius& phaseShiftRadius ) :
    penetrability_( penetrationRadius ),
    shiftFactor_( shiftFactorRadius ),
    phaseShift_( phaseShiftRadius ) {}
