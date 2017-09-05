Type( std::vector< lvalue::Type >&& lvalues,
      double atomicWeightRatio,
      EnergyRange energyRange,
      int nucleonNumber ) :
  lvalues( std::move(lvalues) ),
  target2CompoundWeightRatio( atomicWeightRatio / ( atomicWeightRatio + 1. ) ),
  energyRange( energyRange ),
  nucleonNumber( nucleonNumber ){}
