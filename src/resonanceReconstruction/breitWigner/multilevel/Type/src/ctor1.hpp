Type( std::vector< lvalue::Type >&& lvalues,
      double atomicWeightRatio,
      EnergyRange energyRange,
      int nucleonNumber,
      double targetSpin,
      Radius&& radius ) :
  Parent( std::move(lvalues),
          atomicWeightRatio,
          energyRange,
          nucleonNumber,
          targetSpin ),
  radius( std::move( radius ) ){}
