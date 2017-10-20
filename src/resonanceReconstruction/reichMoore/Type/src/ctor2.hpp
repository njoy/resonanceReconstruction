template< typename Tag >
Type( std::vector< Lvalue >&& lvalues,
      double atomicWeightRatio,
      EnergyRange energyRange,
      Tag tag,
      ChannelRadius&& channelRadius,
      ScatteringRadius&& scatteringRadius ) :
  Parent( std::move(lvalues), atomicWeightRatio, energyRange, tag ),
  channelRadius( std::move( channelRadius ) ),
  scatteringRadius( std::move( scatteringRadius ) ){}
