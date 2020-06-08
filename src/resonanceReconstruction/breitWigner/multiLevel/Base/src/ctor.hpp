Base( std::vector< lvalue::Type >&& lvalues,
      double atomicWeightRatio,
      EnergyRange energyRange,
      int nucleonNumber,
      double targetSpin ) :
  Parent( order( std::move(lvalues) ),
          atomicWeightRatio, energyRange, nucleonNumber ),
  targetSpin( targetSpin ) {}
