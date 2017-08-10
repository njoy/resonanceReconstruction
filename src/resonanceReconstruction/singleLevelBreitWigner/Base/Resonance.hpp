struct Resonance {
  Quantity<ElectronVolts> energy;
  Quantity<ElectronVolts> neutronWidth;
  Quantity<ElectronVolts> captureWidth;
  Quantity<ElectronVolts> fissionWidth;
  Quantity<ElectronVolts> weightedCompetitiveWidth;
  double inversePenetrationFactor;
  double shiftFactor;
  double statisticalFactor;

  CrossSection operator()( Quantity<ElectronVolts> energy,
                           Quantity<InvRootBarns> waveNumber,
                           const double penetrationFactor,
                           const double shiftFactor,
                           const double phaseShift,
                           Quantity<ElectronVolts> competitiveWidth,
                           // from enclosing call
                           const PsiChi& kernel ){      
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
      * this->neutronWidth
      * this->inverseTotalWidth
      / ( waveNumber * waveNumber );

    const auto psichi =
      kernel( energy, primedResonanceEnergy, inverseTotalWidth );
    
    const auto& psi = psichi[0];
    const auto& chi = psichi[1];

    const auto scattering =
      maximumOffset
      * ( ( std::cos( 2. * phaseShift )
            - 1. + neutronWidth * inverseTotalWidth ) * psi
          + std::sin( 2. * phaseShift ) * chi );
      
    const auto capture =
      maximumOffset * this->captureWidth * psi * inverseTotalWidth;

    const auto fission =
      maximumOffset * this->fissionWidth * psi * inverseTotalWidth;
      
    return { scattering, capture, fission };
  }
};
