template< typename PsiChi >
auto operator()( const double penetrationFactor,
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
    this->energy
    + 0.5 * this->neutronWidth
          * this->inversePenetrationFactor
          * ( this->shiftFactor - shiftFactor );
    
  const auto scaling =
    this->statisticalFactor * neutronWidth * inverseTotalWidth;

  const auto psichi =
    kernel( primedResonanceEnergy, inverseTotalWidth );
  
  const auto& psi = psichi[0];
  const auto& chi = psichi[1];

  double scattering =
    scaling * ( ( cos - 1. + neutronWidth * inverseTotalWidth ) * psi + sin * chi );
    
  double capture =
    ( scaling * psi * inverseTotalWidth ) * this->captureWidth;

  double fission =
    ( scaling * psi * inverseTotalWidth ) * this->fissionWidth;
    
  return pack( scattering, capture, fission );
}
