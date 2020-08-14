/**
 *  @brief Constructor
 *
 *  @param[in] energy    the resonance energy (in eV, may be negative)
 *  @param[in] total     the reduced total width
 *  @param[in] elastic   the reduced elastic width
 *  @param[in] capture   the reduced capture width
 *  @param[in] fission   the reduced fission width
 */
Resonance( const Energy& energy,
           const ReducedWidth& total,
           const ReducedWidth& elastic,
           const ReducedWidth& capture,
           const ReducedWidth& fission ) :
    energy_( energy ),
    total_( total ),
    elastic_( elastic ),
    capture_( capture ),
    fission_( fission ),
    delta_( total - capture - fission ) {}
