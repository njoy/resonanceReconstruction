template< typename ChannelRatio >
auto competitiveWidth( const Quantity< ElectronVolts > energy,
                       ChannelRatio&& rho ) const {
  const double penetrationFactor =
    this->penetrationShift( rho( energy + *( this->weightedQValue ) ) )[0];

  return [ penetrationFactor ]( const resonance::Type& r ){
    return r.weightedCompetitiveWidth * penetrationFactor;
  };
}
