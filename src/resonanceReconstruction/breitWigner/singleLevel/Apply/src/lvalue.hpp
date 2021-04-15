template< typename ChannelRatio,
          typename StatisticalFactor,
          typename PenetrationShift >
static auto lvalue( const endf::SingleLevelBreitWigner::LValue& lstate,
                    ChannelRatio&& rho,
                    StatisticalFactor&& g,
                    PenetrationShift&& penetrationShift,
                    double targetSpin ){
  if ( lstate.LRX() ){
    const auto weightRatio = lstate.AWRI() / ( lstate.AWRI() + 1. );
    const auto offset = lstate.QX() * weightRatio * electronVolts;

    auto rs =
      resonances( lstate.resonances(),
                  [ &, GX = competitiveWidth( offset ) ]( auto&& r )
                  { return resonance( r, rho, g, penetrationShift, GX ); } );
    return lvalue::Type( std::move(rs), lstate.L(), offset, targetSpin );
  } else {
    auto rs =
      resonances( lstate.resonances(),
                  [&]( auto&& r )
                  { return resonance( r, rho, g, penetrationShift ); } );
    return lvalue::Type( std::move(rs), lstate.L(), targetSpin );
  }
}

template< typename ChannelRatio,
          typename StatisticalFactor >
static auto lvalue( const endf::SingleLevelBreitWigner::LValue& lstate,
                    ChannelRatio&& rho,
                    StatisticalFactor&& g,
                    double targetSpin ){
  switch( lstate.L() ){
  case 0:
    return lvalue( lstate, rho, g, penetrationShift( Integer<0>{} ), targetSpin );
  case 1:
    return lvalue( lstate, rho, g, penetrationShift( Integer<1>{} ), targetSpin );
  case 2:
    return lvalue( lstate, rho, g, penetrationShift( Integer<2>{} ), targetSpin );
  case 3:
    return lvalue( lstate, rho, g, penetrationShift( Integer<3>{} ), targetSpin );
  case 4:
    return lvalue( lstate, rho, g, penetrationShift( Integer<4>{} ), targetSpin );
  }
  throw std::exception();
}

template< typename BreitWigner,
          typename ChannelRatio,
          typename StatisticalFactor >
static auto lvalues( const BreitWigner& bw,
                     ChannelRatio&& rho,
                     StatisticalFactor&& g,
                     double targetSpin ){
  return
    bw.lValues()
    | ranges::views::transform
    ( [&]( auto&& lstate ){ return lvalue( lstate, rho, g, targetSpin ); } )
    | ranges::to_vector;
}
