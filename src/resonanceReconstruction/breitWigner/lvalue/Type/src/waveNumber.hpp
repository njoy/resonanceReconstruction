auto waveNumber( const Quantity<ElectronVolts> energy ) const {
  using Units4Constant = decltype( pow( Barns() * ElectronVolts(), Ratio<-1,2> ) );

  constexpr Quantity< Units4Constant > neutronConstant =
    5.787793139E-14 * root( kilo(grams) ) / constant::dirac;

  return this->target2CompoundWeightRatio
         * neutronConstant
         * sqrt( std::abs( energy ) );
}
