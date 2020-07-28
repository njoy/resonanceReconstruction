Type( std::vector< lvalue::Type >&& lvalues,
      const double atomicWeightRatio,
      const EnergyRange energyRange,
      const int nucleonNumber ) :
  lvalues( std::move(lvalues) ),
  target2CompoundWeightRatio( atomicWeightRatio / ( atomicWeightRatio + 1. ) ),
  energyRange( energyRange ),
  nucleonNumber( nucleonNumber ){}
