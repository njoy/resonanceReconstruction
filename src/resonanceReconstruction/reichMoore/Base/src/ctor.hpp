template< typename Tag >
Base( std::vector< Lvalue >&& lvalues,
      double atomicWeightRatio,
      EnergyRange energyRange,
      Tag tag ) :
  lvalues( order( std::move(lvalues) ) ),
  energyRange( energyRange ),
  target2CompoundWeightRatio( atomicWeightRatio / ( atomicWeightRatio + 1. ) ),
  tag( tag ){}
