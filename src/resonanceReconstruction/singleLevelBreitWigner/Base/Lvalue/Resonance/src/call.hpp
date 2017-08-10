template< typename PsiChi >
CrossSection operator()( const Quantity<ElectronVolts> energy,
                         const Quantity<InvRootBarns> waveNumber,
                         const double penetrationFactor,
                         const double shiftFactor,
                         const double sin,
                         const double cos,
                         const Quantity<ElectronVolts> competitiveWidth,
                         const PsiChi& kernel ) const {
  const auto neutronWidth =
    this->neutronWidth
    * penetrationFactor
    * this->inversePenetrationFactor;

  const auto inverseTotalWidth =
    1. / ( neutronWidth
           + this->captureWidth
           + this->fissionWidth
           + competitiveWidth );

  const auto primedResonanceEnergy =
    0.5 * this->energy
    + this->neutronWidth
    * this->inversePenetrationFactor
    * ( this->shiftFactor - shiftFactor );
    
  const auto maximumOffset =
    4.0 * pi 
    * this->statisticalFactor
    * neutronWidth
    * inverseTotalWidth
    / ( waveNumber * waveNumber );

  const auto psichi =
    kernel( energy, primedResonanceEnergy, inverseTotalWidth );
  
  const auto& psi = psichi[0];
  const auto& chi = psichi[1];

  const auto scattering =
    maximumOffset
    * ( ( cos - 1. + neutronWidth * inverseTotalWidth ) * psi + sin * chi );
    
  const auto capture =
    maximumOffset * this->captureWidth * psi * inverseTotalWidth;

  const auto fission =
    maximumOffset * this->fissionWidth * psi * inverseTotalWidth;
    
  return { scattering, capture, fission };
}
