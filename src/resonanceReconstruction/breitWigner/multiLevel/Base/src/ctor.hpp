Base( std::vector< lvalue::Type >&& lvalues,
      const double atomicWeightRatio,
      const EnergyRange energyRange,
      const int nucleonNumber,
      const double targetSpin ) :
  Parent( order( std::move(lvalues) ),
          atomicWeightRatio, energyRange, nucleonNumber ),
  targetSpin( targetSpin ) {}
