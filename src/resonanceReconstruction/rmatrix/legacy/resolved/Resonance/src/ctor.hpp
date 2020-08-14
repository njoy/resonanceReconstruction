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
           const Width& total,
           const Width& elastic,
           const Width& capture,
           const Width& fission,
           double penetrability,
           double shiftfactor ) :
    energy_( energy ),
    total_( total ),
    elastic_( elastic ),
    capture_( capture ),
    fission_( fission ),
    elastic_to_penetrability_( elastic / penetrability ),
    penetrability_( penetrability ),
    shiftfactor_( shiftfactor ),
    delta_( total - capture - fission ) {}
