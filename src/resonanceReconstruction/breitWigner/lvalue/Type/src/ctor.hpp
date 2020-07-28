Type( std::vector< resonance::Type >&& resonances,
      const int orbitalAngularMomentum,
      const double targetSpin ) :
  resonances( std::move( resonances ) ),
  weightedQValue( std::nullopt ),
  orbitalAngularMomentum( orbitalAngularMomentum ),
  statisticalFactorSum( statisticalFactorSummation
                        ( targetSpin, orbitalAngularMomentum ) ){}

Type( std::vector< resonance::Type >&& resonances,
      const int orbitalAngularMomentum,
      const Quantity< ElectronVolts > weightedQValue,
      const double targetSpin ) :
  resonances( std::move( resonances ) ),
  weightedQValue( weightedQValue ),
  orbitalAngularMomentum( orbitalAngularMomentum ),
  statisticalFactorSum( statisticalFactorSummation
                        ( targetSpin, orbitalAngularMomentum ) ){}
