template< typename PsiChi >
auto operator()( const double penetrationFactor,
                 const double shiftFactor,
                 const double sin,
                 const double cos,
                 const PsiChi& kernel,
                 const Quantity<ElectronVolts>
                 const competitiveWidth = 0.0 * barns ) const {

  const auto weightedNeutronWidth =
    this->neutronWidth * this->inversePenetrationFactor;

  const auto neutronWidth = weightedNeutronWidth * penetrationFactor;

  const auto inverseTotalWidth =
    1. / ( neutronWidth
           + this->captureWidth
           + this->fissionWidth
           + competitiveWidth );

  const auto primedResonanceEnergy =
    this->energy
    + 0.5 * weightedNeutronWidth * ( this->shiftFactor - shiftFactor );

  const auto psichi = kernel( primedResonanceEnergy, inverseTotalWidth );

  const auto& psi = psichi[0];
  const auto& chi = psichi[1];

  const double widthRatio = neutronWidth * inverseTotalWidth;
  const double scattering = ( ( cos - 1. + widthRatio ) * psi + sin * chi );
  const double factor = psi * inverseTotalWidth;
  const double capture = factor * this->captureWidth;
  const double fission = factor * this->fissionWidth;

  return
    pack( scattering, capture, fission ) * this->statisticalFactor * widthRatio;
}
