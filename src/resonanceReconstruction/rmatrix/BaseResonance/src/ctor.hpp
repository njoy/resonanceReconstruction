/**
 *  @brief Constructor
 *
 *  @param[in] energy       the resonance energy (in eV, may be negative)
 *  @param[in] widths       the reduced widths (in sqrt(eV))
 */
BaseResonance( const Energy& energy,
               std::vector< ReducedWidth >&& widths ) :
    energy_( energy ),
    widths_( std::move( widths ) ) {}
