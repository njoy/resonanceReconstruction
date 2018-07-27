/**
 *  @brief Constructor
 *
 *  @param[in] energy       the resonance energy (in eV, may be negative)
 *  @param[in] widths       the reduced widths (in sqrt(eV))
 *  @param[in] eliminated   the reduced eliminated capture width (in sqrt(eV))
 */
Resonance( const Energy& energy,
           std::vector< ReducedWidth >&& widths,
           const ReducedWidth& eliminated = 0.0 * rootElectronVolt ) :
    energy_( energy ),
    widths_( std::move( widths ) ),
    eliminated_( eliminated ) {}

