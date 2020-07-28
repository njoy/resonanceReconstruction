Type( std::vector< lvalue::Type >&& lvalues,
      const double atomicWeightRatio,
      const EnergyRange energyRange,
      const int nucleonNumber,
      ChannelRadius&& channelRadius,
      ScatteringRadius&& scatteringRadius ) :
  Parent( std::move(lvalues),
          atomicWeightRatio,
          energyRange,
          nucleonNumber ),
  channelRadius( std::move( channelRadius ) ),
  scatteringRadius( std::move( scatteringRadius ) ){}
