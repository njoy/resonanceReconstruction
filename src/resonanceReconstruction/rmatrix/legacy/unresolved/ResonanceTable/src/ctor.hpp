private:

/**
 *  @brief Private intermediate constructor
 */
ResonanceTable( std::vector< Resonance >&& resonances,
                Degrees&& degrees,
                std::tuple< LevelSpacingTable, ReducedWidthTable,
                            WidthTable, WidthTable, WidthTable >&& tables ) :
    ResonanceTableBase( std::move( resonances ) ),
    degrees_( std::move( degrees ) ),
    level_spacing_table_( std::move( std::get< 0 >( tables ) ) ),
    elastic_table_( std::move( std::get< 1 >( tables ) ) ),
    capture_table_( std::move( std::get< 2 >( tables ) ) ),
    fission_table_( std::move( std::get< 3 >( tables ) ) ),
    competition_table_( std::move( std::get< 4 >( tables ) ) ) {

    verifyTable( this->resonances(), this->degrees_ );
}

public:

/**
 *  @brief Constructor
 *
 *  @param[in] resonances   the resonances (ne values)
 *  @param[in] widths       the degrees of freedom for each channel
 *                          (nc values between 0 an 4)
 */
ResonanceTable( std::vector< Resonance >&& resonances,
                Degrees&& degrees ) :
    ResonanceTable( std::move( resonances ),
                    std::move( degrees ),
                    makeTables( resonances ) ) {}
