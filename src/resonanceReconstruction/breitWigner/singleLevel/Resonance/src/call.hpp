template< typename PsiChi >
auto operator()( const double penetrationFactor,
                 const double shiftFactor,
                 const double sinSquaredPhi,
                 const double sin2Phi,
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

  const double scattering =
    ( widthRatio - 2 * sinSquaredPhi ) * psi + sin2Phi * chi;

  const auto scaling = psi * inverseTotalWidth;
  const double capture = scaling * this->captureWidth;
  const double fission = scaling * this->fissionWidth;

  return pack( scattering, capture, fission )
         * ( this->statisticalFactor * widthRatio );
}
