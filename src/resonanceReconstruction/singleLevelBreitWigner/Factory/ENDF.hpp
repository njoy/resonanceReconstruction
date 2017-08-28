template< >
class Factory< ENDFtk::resonanceParameters::resolved::SLBW > {
  namespace ENDF = ENDFtk::resonanceParameters::resolved;
  
  template< typename ChannelRadius >
  static auto lvalue( const ENDF::SLBW& slbw,
                      const ENDFtk::LIST& list,
                      ChannelRadius&& radius ){

    const double target2CompoundWeightRatio =
      list.C1() / ( list.C1() + 1. );

    const std::optional< Quantity<ElectronVolts > >
      weightedQvalue = ( slbw.LRX() ) ?
      list.C2() * electronVolts / target2CompoundWeightRatio :
      std::nullopt;
    
    const auto processResonance = [&]( auto&& chunk ){
      const auto energy = chunk[0] * electronVolts;
      const auto neutronWidth = chunk[3] * electronVolts;
      const auto captureWidth = chunk[4] * electronVolts;
      const auto fissionWidth = chunk[5] * electronVolts;
      
      const auto waveNumber = neutronWaveNumber( list.C1() );
  
      const auto weightedCompetitiveWidth =
        ( not slbw.LRX() ) ? 0.0 * electronVolts
                           : [energy = energy + *weightedQvalue, & ]{
        const auto channelRadius = radius( std::abs(energy) );
        const auto channelRatio = channelRadius * waveNumber( energy );
      
        const auto penetrationFactor =
          penetrationShift( list.L1(), channelRatio )[0];
        
        return ranges::accumulate( chuck | ranges::view::slice(3,6),
                                   chunk[2],
                                   ranges::minus{} ) / penetrationFactor;
      }();
  
      const auto channelRadius = radius( std::abs(energy) );
      const auto channelRatio = channelRadius * waveNumber( energy );
      
      const auto ps = penetrationShift( list.L1(), channelRatio )[0];
      const auto& penetrationFactor = ps[0];
      const auto& shiftFactor = ps[1];
      const auto statisticalFactor = ( 2 * chunk[1] + 1 )
                                     / ( 4 * slbw.SPI() + 2 );
      return Base::Lvalue::Resonance{ energy,
                                      neutronWidth,
                                      captureWidth,
                                      fissionWidth,
                                      weightedCompetitiveWidth,
                                      1. / penetrationFactor,
                                      shiftFactor,
                                      statisticalFactor };
    };

    auto resonances = list.B()
                      | ranges::view::chunk(6)
                      | ranges::view::transform( processResonance )
                      | ranges::to_vector;

    return Base::Lvalue{ std::move(resonances),
                         weightedQvalue,
                         target2CompoundWeightRatio,
                         list.L1() };
  }
  
  template< typename Lists >
  static auto lvalues( const ENDF::SLBW& slbw,
                       ChannelRadius&& radius ){
    return slbw.LISTS()
           | ranges::view::transform
             ( [&]( auto&& list ){ return lvalue( slbw, list, radius ); })
           | ranges::to_vector;
  }

public: 
  
  template< typename Radius >
  static auto build( const ENDF::SLBW& slbw, Radius&& radius ){
    auto ls = lvalues( slbw, radius );
    return Type<Radius>{ std::move(radius),
                         std::move(ls),
                         Base::EnergyRange{ slbw.EL() * electronVolts,
                                            slbw.EH() * electronVolts },
                         std::ceil( slbw.LISTS().front().C1() ) };
  }

  template< typename ChannelRadius, typename ScatteringRadius >
  static auto build( const ENDF::SLBW& slbw,
                     ChannelRadius&& channelRadius,
                     ScatteringRadius&& scatteringRadius ){
    auto ls = lvalues( slbw, channelRadius );
    return Type< ChannelRadius, ScatteringRadius >
           { std::move(channelRadius),
             std::move(scatteringRadius),
             std::move(ls),
             Base::EnergyRange{ slbw.EL() * electronVolts,
                                slbw.EH() * electronVolts },
               std::ceil( slbw.LISTS().front().C1() ) };
  }

  template< typename Functor >
  decltype(auto) apply( const ENDF::SLBW& slbw, Functor&& callback ) const {
    switch( slbw.NRO() ){
    case 0:
      switch( slbw.NAPS() ){
      case 0: {
        return callback( build( slbw,
                                channelRadius( slbw.LISTS().front().C1() ),
                                radius( slbw.AP() ) ) );
      }
      case 1:{
        return callback( build( slbw,
                                radius( slbw.AP() ) ) );
      }
      }
    case 1:
      switch( slbw.NAPS() ){
      case 0:{
        return callback( build( slbw,
                                channelRadius( slbw.LISTS().front().C1() ),
                                radius( slbw.APE() ) ) );
      }
      case 1:{
        return callback( build( slbw,
                                radius( slbw.APE() ) ) );
      }
      case 2:{
        return callback( build( slbw,
                                radius( slbw.AP() ),
                                radius( slbw.APE() ) ) );
      }
      }
    } 
  }
};
