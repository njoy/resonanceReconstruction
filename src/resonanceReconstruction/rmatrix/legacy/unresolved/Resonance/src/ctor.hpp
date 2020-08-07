/**
 *  @brief Constructor
 *
 *  @param[in] energy       the resonance energy (in eV, may be negative)
 *  @param[in] eliminated   the level spacing energy (in eV)
 *  @param[in] widths       the reduced widths (in sqrt(eV))
 */
Resonance( const Energy& energy,
           const Energy& spacing,
           const ReducedWidth& elastic,
           const Width& capture,
           const Width& fission,
           const Width& competition ) :
    energy_( energy ),
    spacing_( spacing ),
    elastic_( elastic ),
    capture_( capture ),
    fission_( fission ),
    competition_( competition ) {}
