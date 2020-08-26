/**
 *  @brief Constructor
 *
 *  @param[in] table   the energy dependent channel radius
 */
 ChannelRadiusTable( RadiusTable&& table ) :
   table_( std::make_shared< RadiusTable >( std::move( table ) ) ) {}

/**
 *  @brief Default copy constructor
 *
 *  @param[in] table   the instance to be copied
 */
ChannelRadiusTable( const ChannelRadiusTable& ) = default;

/**
 *  @brief Default copy assignment
 *
 *  @param[in] table   the instance to be copied
 */
ChannelRadiusTable& operator=( const ChannelRadiusTable& ) = default;
