template< typename Tag, typename Radius >
static auto
build( const ENDF::resolved::ReichMoore& rm,
       Tag tag,
       Radius&& radius,
       bool useAPL ){
  const double atomicWeightRatio = rm.lValues().front().AWRI();

  const auto k = neutronWaveNumber( atomicWeightRatio );

  const auto g =
    [ denominator = 1. / ( 4 * rm.SPI() + 2. ) ]
    ( const double J ){ return ( 2 * J + 1. ) * denominator; };

  auto lstates = lvalues( rm, k, radius, g, useAPL );

  return Type< std::decay_t< Radius > >
    ( std::move(lstates),
      atomicWeightRatio,
      EnergyRange{ rm.EL() * electronVolts,
                   rm.EH() * electronVolts },
      tag,
      std::move(radius) );
}

template< typename Tag, typename ChannelRadius, typename ScatteringRadius >
static auto build( const ENDF::resolved::ReichMoore& rm,
                   Tag tag,
                   ChannelRadius&& channelRadius,
                   ScatteringRadius&& scatteringRadius,
                   bool useAPL ){

  const double atomicWeightRatio = rm.lValues().front().AWRI();

  const auto k = neutronWaveNumber( atomicWeightRatio );

  const auto g =
    [ denominator = 1. / ( 4 * rm.SPI() + 2. ) ]
    ( const double J ){ return ( 2 * J + 1. ) * denominator; };

  auto lstates = lvalues( rm, k, channelRadius, g, useAPL );

  return Type< std::decay_t< ChannelRadius >,
               std::decay_t< ScatteringRadius > >
    ( std::move(lstates),
      atomicWeightRatio,
      EnergyRange{ rm.EL() * electronVolts,
                   rm.EH() * electronVolts },
      tag,
      std::move(channelRadius),
      std::move(scatteringRadius) );
}
