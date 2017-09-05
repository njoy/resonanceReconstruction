Type( std::vector< lvalue::Type >&& lvalues,
      double atomicWeightRatio,
      EnergyRange energyRange,
      int nucleonNumber,
      Radius&& radius ) :
  Base( std::move(lvalues),
        atomicWeightRatio,
        energyRange,
        nucleonNumber ),
  radius( std::move( radius ) ){}
