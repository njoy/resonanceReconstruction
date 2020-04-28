template< typename Radius >
static auto build( const EnergyRange& energyRange,
                   const ENDF::resolved::SLBW& slbw,
                   Radius&& radius ){
  const double atomicWeightRatio = slbw.lValues().front().AWRI();

  const auto rho =
    [ &, k = neutronWaveNumber( atomicWeightRatio ) ]
    ( const Quantity< ElectronVolts > energy ){
    return k( energy ) * radius( energy );
  };

  const auto g =
    [ denominator = 1. / ( 4 * slbw.SPI() + 2. ) ]
    ( const double J ){ return ( 2 * J + 1. ) * denominator; };

  auto lstates = lvalues( slbw, rho, g, slbw.SPI() );

  return Type< std::decay_t< Radius > >
    ( std::move(lstates),
      atomicWeightRatio,
      energyRange,
      nucleonNumber( slbw ),
      std::move(radius) );
}

template< typename ChannelRadius, typename ScatteringRadius >
static auto build( const EnergyRange& energyRange,
                   const ENDF::resolved::SLBW& slbw,
                   ChannelRadius&& channelRadius,
                   ScatteringRadius&& scatteringRadius ){

  const double atomicWeightRatio = slbw.lValues().front().AWRI();

  const auto rho =
    [ &, k = neutronWaveNumber( atomicWeightRatio ) ]
    ( const Quantity< ElectronVolts > energy ){
    return k( energy ) * channelRadius( energy );
  };

  const auto g =
    [ denominator = 1. / ( 4 * slbw.SPI() + 2. ) ]
    ( const double J ){ return ( 2 * J + 1. ) * denominator; };

  auto lstates = lvalues( slbw, rho, g, slbw.SPI() );

  return Type< std::decay_t< ChannelRadius >,
               std::decay_t< ScatteringRadius > >
    ( std::move(lstates),
      atomicWeightRatio,
      energyRange,
      nucleonNumber( slbw ),
      std::move(channelRadius),
      std::move(scatteringRadius) );
}
