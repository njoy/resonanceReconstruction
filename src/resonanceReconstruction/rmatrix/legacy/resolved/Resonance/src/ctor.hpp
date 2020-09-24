/**
 *  @brief Constructor
 *
 *  @param[in] energy    the resonance energy Er (in eV, may be negative)
 *  @param[in] elastic   the reduced elastic width
 *  @param[in] capture   the reduced capture width
 *  @param[in] fission   the reduced fission width
 *  @param[in] total     the competitive total width
 *  @param[in] P         the penetrability at the energy Er
 *  @param[in] Q         the penetrability at the energy Er - Q
 *  @param[in] S         the shift factor at the energy Er
 */
Resonance( const Energy& energy,
           const Width& elastic,
           const Width& capture,
           const Width& fission,
           const Width& competition,
           double P,
           double Q,
           double S ) :
    energy_( energy ),
    elastic_( elastic ),
    capture_( capture ),
    fission_( fission ),
    competition_( competition ),
    elastic_to_penetrability_( elastic / P ),
    competition_to_penetrability_( competition / Q ),
    shiftfactor_( S ) {}
