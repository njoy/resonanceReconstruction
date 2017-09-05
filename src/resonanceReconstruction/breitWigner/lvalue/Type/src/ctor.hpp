Type( std::vector< resonance::Type >&& resonances,
      int orbitalAngularMomentum ) :
  resonances( std::move( resonances ) ),
  weightedQValue( std::nullopt ),
  orbitalAngularMomentum( orbitalAngularMomentum ){}

Type( std::vector< resonance::Type >&& resonances,
      int orbitalAngularMomentum,
      Quantity<ElectronVolts > weightedQValue ) :
  resonances( std::move( resonances ) ),
  weightedQValue( weightedQValue ),
  orbitalAngularMomentum( orbitalAngularMomentum ){} 
