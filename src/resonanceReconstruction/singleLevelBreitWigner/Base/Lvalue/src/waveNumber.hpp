auto waveNumber( Quantity<ElectronVolts> energy ) const {
  constexpr auto neutronMass = 1.674927471E-21 * kilo(grams);

  using Units4Constant =
    decltype( pow( Barns() * ElectronVolts(), Ratio<-1,2> ) );
  
  const Quantity< Units4Constant >
    neutronConstant = sqrt( 2.0 * neutronMass ) / constant::dirac;
  
  return this->target2CompoundWeightRatio
         * neutronConstant
         * sqrt( std::abs( energy ) );
}
