/**
 *  @brief Constructor
 *
 *  @param[in] energy       the resonance energy (in eV, may be negative)
 *  @param[in] eliminated   the level spacing energy (in eV)
 *  @param[in] widths       the reduced widths (in sqrt(eV))
 */
UnresolvedResonance( const Energy& energy,
                     const Energy& spacing,
                     std::vector< ReducedWidth >&& widths ) :
    BaseResonance( energy, std::move( widths ) ),
    spacing_( spacing ) {}
