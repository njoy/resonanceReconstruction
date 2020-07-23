/**
 *  @brief Constructor
 *
 *  Both energy values and reduced widths may be negative. When the energy is
 *  negative, the penetrability, shift factor, phase shift, etc. are calculated
 *  at the absolute value of the energy.
 *
 *  When using this Resonance class with general R matrix theory, the reduced
 *  eliminated capture width does not need to be specified (it has a default
 *  value of zero).
 *
 *  @param[in] energy       the resonance energy (in eV, may be negative)
 *  @param[in] widths       the reduced widths (in sqrt(eV))
 *  @param[in] eliminated   the reduced eliminated capture width (in sqrt(eV))
 */
Resonance( const Energy& energy,
           std::vector< ReducedWidth >&& widths,
           const ReducedWidth& eliminated = 0.0 * rootElectronVolt ) :
    BaseResonance( energy, std::move( widths ) ),
    eliminated_( eliminated ) {}
