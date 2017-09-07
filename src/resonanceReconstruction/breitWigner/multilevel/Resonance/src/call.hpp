template< typename PsiChi >
auto operator()( const double penetrationFactor,
                 const double shiftFactor,
                 const Quantity< ElectronVolts > competitiveWidth,
                 const PsiChi& kernel ) const {
  const auto weightedWidth =
    this->neutronWidth * this->inversePenetrationFactor;
  
  const auto neutronWidth = weightedWidth * penetrationFactor;

  const auto inverseTotalWidth =
    1. / ( neutronWidth
           + this->captureWidth
           + this->fissionWidth
           + competitiveWidth );

  const auto primedResonanceEnergy =
    this->energy
    + 0.5 * weightedWidth * ( this->shiftFactor - shiftFactor );

  const auto psichi = kernel( primedResonanceEnergy, inverseTotalWidth );
  
  const double& psi = psichi[0];
  const double& chi = psichi[1];
  
  const double widthRatio = neutronWidth * inverseTotalWidth;
  const double scattering = pack( psi, chi ) * widthRatio;

  const auto scaling = psi * inverseTotalWidth * widthRatio;
  const double capture = scaling * this->captureWidth;
  const double fission = scaling * this->fissionWidth;
    
  return pack( scattering, capture, fission ); 
}
