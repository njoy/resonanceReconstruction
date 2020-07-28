Type( std::vector< lvalue::Type >&& lvalues,
      const double atomicWeightRatio,
      const EnergyRange energyRange,
      const int nucleonNumber,
      Radius&& radius ) :
  Parent( std::move(lvalues),
          atomicWeightRatio,
          energyRange,
          nucleonNumber ),
  radius( std::move( radius ) ){}
