template< typename Radius >
static auto build( const ENDF::resolved::MLBW& mlbw, Radius&& radius ){
  const double atomicWeightRatio = mlbw.lValues().front().AWRI();

  const auto rho =
    [ &, k = neutronWaveNumber( atomicWeightRatio ) ]
    ( const Quantity< ElectronVolts > energy ){
    return k( energy ) * radius( energy );
  };

  const auto g =
    [ denominator = 1. / ( 4 * mlbw.SPI() + 2. ) ]
    ( const double J ){ return ( 2 * J + 1. ) * denominator; };

  auto lstates = lvalues( mlbw, rho, g, mlbw.SPI() );

  return Type< std::decay_t< Radius > >
    ( std::move(lstates),
      atomicWeightRatio,
      EnergyRange{ mlbw.EL() * electronVolts,
                   mlbw.EH() * electronVolts },
      nucleonNumber( mlbw ),
      mlbw.SPI(),
      std::move(radius) );
}

template< typename ChannelRadius, typename ScatteringRadius >
static auto build( const ENDF::resolved::MLBW& mlbw,
                   ChannelRadius&& channelRadius,
                   ScatteringRadius&& scatteringRadius ){

  const double atomicWeightRatio = mlbw.lValues().front().AWRI();

  const auto rho =
    [ &, k = neutronWaveNumber( atomicWeightRatio ) ]
    ( const Quantity< ElectronVolts > energy ){
    return k( energy ) * channelRadius( energy );
  };

  const auto g =
    [ denominator = 1. / ( 4 * mlbw.SPI() + 2. ) ]
    ( const double J ){ return ( 2 * J + 1. ) * denominator; };

  auto lstates = lvalues( mlbw, rho, g, mlbw.SPI() );

  return Type< std::decay_t< ChannelRadius >,
               std::decay_t< ScatteringRadius > >
    ( std::move(lstates),
      atomicWeightRatio,
      EnergyRange{ mlbw.EL() * electronVolts,
                   mlbw.EH() * electronVolts },
      nucleonNumber( mlbw ),
      mlbw.SPI(),
      std::move(channelRadius),
      std::move(scatteringRadius) );
}
