Type( std::vector< lvalue::Type >&& lvalues,
      const double atomicWeightRatio,
      const EnergyRange energyRange,
      const int nucleonNumber,
      const double targetSpin,
      ChannelRadius&& channelRadius,
      ScatteringRadius&& scatteringRadius ) :
  Parent( std::move(lvalues),
          atomicWeightRatio,
          energyRange,
          nucleonNumber,
          targetSpin ),
  channelRadius( std::move( channelRadius ) ),
  scatteringRadius( std::move( scatteringRadius ) ){}
