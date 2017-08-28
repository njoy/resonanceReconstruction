template< typename T >
constexpr auto root( T t ){
  using Root = decltype( pow( typename decltype( 1.0 * t )::Units(), Ratio<1,2> ) );
  Quantity< Root > r; r.value = 1.0;
  return r;
}

inline auto neutronWaveNumber( const double atomicWeightRatio ){
  /* sqrt( 2.0 * neutron mass ) */
  using Units4Constant =
    decltype( pow( Barns() * ElectronVolts(), Ratio<-1,2> ) );
  
  constexpr Quantity< Units4Constant >
    neutronConstant = 5.787793139E-14 * root( kilo(grams) ) / constant::dirac;
  
  const auto weightFraction =
    atomicWeightRatio / ( atomicWeightRatio + 1.0 );

  return
    [ constant = weightFraction * neutronConstant ]
    ( Quantity<ElectronVolts> energy )
    {
      return constant * sqrt( std::abs( energy ) );
    };
}
