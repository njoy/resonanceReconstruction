Type( std::vector< lvalue::Type >&& lvalues,
      double atomicWeightRatio,
      EnergyRange energyRange,
      int nucleonNumber,
      ChannelRadius&& channelRadius,
      ScatteringRadius&& scatteringRadius ) :
  Parent( std::move(lvalues),
          atomicWeightRatio,
          energyRange,
          nucleonNumber ),
  channelRadius( std::move( channelRadius ) ),
  scatteringRadius( std::move( scatteringRadius ) ){}
