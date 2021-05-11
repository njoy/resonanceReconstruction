//! @todo pybind11 variant needs default constructor workaround
#ifdef PYBIND11
/**
 *  @brief Default constructor - only enabled for pybind11
 */
ChannelRadii() = default;
#endif

/**
 *  @brief Constructor
 *
 *  @param[in] radius   the channel radius to be used for P, S and phi
 */
ChannelRadii( const ChannelRadiusVariant& radius ) :
    penetrability_( radius ),
    shiftFactor_( radius ),
    phaseShift_( radius ) {}

/**
 *  @brief Constructor
 *
 *  @param[in] trueRadius        the channel radius to be used for P and S
 *  @param[in] effectiveRadius   the channel radius to be used for phi
 */
ChannelRadii( const ChannelRadiusVariant& trueRadius,
              const ChannelRadiusVariant& effectiveRadius ) :
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
ChannelRadii( const ChannelRadiusVariant& penetrationRadius,
              const ChannelRadiusVariant& shiftFactorRadius,
              const ChannelRadiusVariant& phaseShiftRadius ) :
    penetrability_( penetrationRadius ),
    shiftFactor_( shiftFactorRadius ),
    phaseShift_( phaseShiftRadius ) {}
