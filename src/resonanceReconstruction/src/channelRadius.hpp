inline auto channelRadius( double atomicWeightRatio ){
  constexpr double amu = 1.00866491582;
  const auto scalar =
    0.123 * std::pow( amu * atomicWeightRatio, 1. / 3. ) + 0.08;
  
  return
    [ radius = scalar * rootBarn ]
    ( Quantity< ElectronVolts > /* energy */ )
    { return radius; };
}
