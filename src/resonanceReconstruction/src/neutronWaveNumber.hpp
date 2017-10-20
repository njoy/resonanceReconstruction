inline auto neutronWaveNumber( const double atomicWeightRatio ){
  using Units4Constant =
    decltype( pow( Barns() * ElectronVolts(), Ratio<-1,2> ) );

  constexpr Quantity< Units4Constant >
                      /* sqrt( 2.0 * neutron mass ) */
    neutronConstant = 5.787793139E-14 * root( kilo(grams) ) / constant::dirac;

  const auto weightFraction =
    atomicWeightRatio / ( atomicWeightRatio + 1.0 );

  return
    [ constant = weightFraction * neutronConstant ]
    ( const Quantity<ElectronVolts> energy )
    { return constant * sqrt( std::abs( energy ) ); };
}
