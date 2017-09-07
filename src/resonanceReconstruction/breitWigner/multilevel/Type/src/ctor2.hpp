Type( std::vector< lvalue::Type >&& lvalues,
      double atomicWeightRatio,
      EnergyRange energyRange,
      int nucleonNumber,
      double targetSpin,
      ChannelRadius&& channelRadius,
      ScatteringRadius&& scatteringRadius ) :
  Base( std::move(lvalues),
        atomicWeightRatio,
        energyRange,
        nucleonNumber,
        targetSpin ),
  channelRadius( std::move( channelRadius ) ),
  scatteringRadius( std::move( scatteringRadius ) ){}
