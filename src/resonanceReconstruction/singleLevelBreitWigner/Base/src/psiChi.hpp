static auto
psiChi( const Quantity<ElectronVolts> energy ){
  return [ energy ]( const Quantity<ElectronVolts> primedResonanceEnergy,
                     const Quantity<InvElectronVolts> inverseTotalWidth ) {
    const double x =
      2. * ( energy - primedResonanceEnergy ) * inverseTotalWidth;
    
    const auto psi = 1. / ( 1. + ( x * x ) );
    return std::array< double, 2 >{{ psi, x * psi }};
  };
}
