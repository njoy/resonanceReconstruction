Type( std::vector< lvalue::Type >&& lvalues,
      double atomicWeightRatio,
      EnergyRange energyRange,
      int nucleonNumber,
      double targetSpin,
      Radius&& radius ) :
  Base( std::move(lvalues),
        atomicWeightRatio,
        energyRange,
        nucleonNumber,
        targetSpin ),
  radius( std::move( radius ) ){}
