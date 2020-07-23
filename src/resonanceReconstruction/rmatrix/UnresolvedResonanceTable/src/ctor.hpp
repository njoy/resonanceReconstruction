private:

/**
 *  @brief Private intermediate constructor
 */
UnresolvedResonanceTable( std::vector< ChannelID >&& channels,
                          std::vector< UnresolvedResonance >&& resonances,
                          std::vector< unsigned int >&& degrees,
                          Table< Energy, Energy > spacingTable,
                          std::vector< Table< Energy, ReducedWidth > > tables ) :
    BaseResonanceTable( std::move( channels ), std::move( resonances ) ),
    degrees_( std::move( degrees ) ),
    level_spacing_table_( std::move( spacingTable ) ),
    width_tables_( std::move( tables ) ) {

    verifyTable( this->channels(), this->degrees_ );
}

public:

/**
 *  @brief Constructor
 *
 *  @param[in] channels     the channel IDs (nc values)
 *  @param[in] resonances   the resonances (ne values)
 *  @param[in] widths       the degrees of freedom for each channel
 *                          (nc values between 0 an 4)
 */
UnresolvedResonanceTable( std::vector< ChannelID >&& channels,
                          std::vector< UnresolvedResonance >&& resonances,
                          std::vector< unsigned int >&& degrees ) :
    UnresolvedResonanceTable( std::move( channels ), std::move( resonances ),
                              std::move( degrees ),
                              makeLevelSpacingTable( resonances ),
                              makeWidthTables( resonances ) ) {}
