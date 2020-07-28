Type( std::vector< lvalue::Type >&& lvalues,
      const double atomicWeightRatio,
      const EnergyRange energyRange,
      const int nucleonNumber,
      const double targetSpin,
      Radius&& radius ) :
  Parent( std::move(lvalues),
          atomicWeightRatio,
          energyRange,
          nucleonNumber,
          targetSpin ),
  radius( std::move( radius ) ){}
