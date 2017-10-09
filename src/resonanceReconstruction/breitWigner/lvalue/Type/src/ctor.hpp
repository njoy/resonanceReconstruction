Type( std::vector< resonance::Type >&& resonances,
      int orbitalAngularMomentum,
      double targetSpin ) :
  resonances( std::move( resonances ) ),
  weightedQValue( std::nullopt ),
  orbitalAngularMomentum( orbitalAngularMomentum ),
  statisticalFactorSum( statisticalFactorSummation
                        ( targetSpin, orbitalAngularMomentum ) ){}

Type( std::vector< resonance::Type >&& resonances,
      int orbitalAngularMomentum,
      Quantity< ElectronVolts > weightedQValue,
      double targetSpin ) :
  resonances( std::move( resonances ) ),
  weightedQValue( weightedQValue ),
  orbitalAngularMomentum( orbitalAngularMomentum ),
  statisticalFactorSum( statisticalFactorSummation
                        ( targetSpin, orbitalAngularMomentum ) ){}
