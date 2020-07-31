Quantity< InvRootBarns >
neutronWaveNumber( const Quantity< ElectronVolts > energy ) const {
  using Units4Constant =
    decltype( pow( Barns() * ElectronVolts(), Ratio<-1,2> ) );

  constexpr Quantity< Units4Constant >
    // sqrt( 2. * neutronMass ) / constant::dirac;
    neutronConstant = 5.787793139E-14 * root( kilo(grams) ) / constant::dirac;

  return this->target2CompoundWeightRatio
         * neutronConstant
         * sqrt( std::abs( energy ) );
}

auto neutronWaveNumber() const {
   return [this]( const Quantity< ElectronVolts > energy )
   { return neutronWaveNumber(energy); };
}
