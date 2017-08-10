auto operator()( const Quantity<ElectronVolts>& energy,
                 const Quantity<InvRootBarns>& waveNumber,
                 const double& penetrationFactor,
                 const double& shiftFactor,
                 const double& phaseShift,
                 const Quantity<ElectronVolts>& competitiveWidth,
                 const PsiChi& kernel ) const {

  const auto process = [&]( auto&& rEnergy,
                            auto&& rNeutronWidth,
                            auto&& rCaptureWidth,
                            auto&& rFissionWidth,
                            auto&& rInversePenetrationFactor,
                            auto&& rShiftFactor,
                            auto&& rStatisticalFactor )-> CrossSection {
    const auto neutronWidth =
      penetrationFactor
      * rNeutronWidth
      * rInversePenetrationFactor;

    const auto inverseTotalWidth =
      1. / ( neutronWidth
             + competitiveWidth
             + rCaptureWidth
             + rFissionWidth );

    const auto primedResonanceEnergy =
      0.5 * rEnergy
      + rNeutronWidth
        * rInversePenetrationFactor
        * ( rShiftFactor - shiftFactor );
  
    const auto maximumOffset =
      4.0 * pi 
      * rStatisticalFactor
      * neutronWidth
      * inverseTotalWidth
      / ( waveNumber * waveNumber );

    const auto psichi = kernel( energy,
                                primedResonanceEnergy,
                                inverseTotalWidth );

    const auto& psi = psichi[0];
    const auto& chi = psichi[1];

    const auto scattering =
      maximumOffset
      * ( ( std::cos( 2. * phaseShift )
            - 1. + neutronWidth * inverseTotalWidth ) * psi
          + std::sin( 2. * phaseShift ) * chi );
  
    const auto capture =
      maximumOffset * rCaptureWidth * psi * inverseTotalWidth;

    const auto fission =
      maximumOffset * rFissionWidth * psi * inverseTotalWidth;

    return { scattering, capture, fission };
  }
  
  return ranges::view::zip_with( process,
                                 this->energy(),
                                 this->neutronWidth(),
                                 this->captureWidth(),
                                 this->fissionWidth(),
                                 this->inversePenetrationFactor(),
                                 this->shiftFactor(),
                                 this->statisticalFactor() );
}
