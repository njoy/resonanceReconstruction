template< typename Tag >
Type( std::vector< Lvalue >&& lvalues,
      double atomicWeightRatio,
      EnergyRange energyRange,
      Tag tag,
      Radius&& radius ) :
  Parent( std::move(lvalues), atomicWeightRatio, energyRange, tag ),
  radius( std::move( radius ) ){}
