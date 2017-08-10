inline auto neutronWaveNumber( double atomicWeightRatio ){
  constexpr auto neutronMass = 1.674927471E-21 * kilo(grams);

  using Units4Constant =
    decltype( pow( Barns() * ElectronVolts(), Ratio<-1,2> ) );
  
  const Quantity< Units4Constant >
    neutronConstant = sqrt( 2.0 * neutronMass ) / constant::dirac;
  
  const auto weightFraction =
    atomicWeightRatio / ( atomicWeightRatio + 1.0 );

  return
    [ constant = weightFraction * neutronConstant ]
    ( Quantity<ElectronVolts> energy )
    { return constant * sqrt( std::abs( energy ) ); };
}
