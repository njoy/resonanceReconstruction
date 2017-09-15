template< typename ChannelRatio,
          typename StatisticalFactor,
          typename PenetrationShift >
static auto lvalue( const ENDF::resolved::SLBW::LState& lstate,
                    ChannelRatio&& rho,
                    StatisticalFactor&& g,
                    PenetrationShift&& penetrationShift ){
  if ( lstate.LRX() ){
    const auto weightRatio = lstate.AWRI() / ( lstate.AWRI() + 1. );
    const auto offset = lstate.QX() * weightRatio * electronVolts;
    
    auto rs =
      resonances( lstate.resonances(),
                  [ &, GX = competitiveWidth( offset ) ]( auto&& r )
                  { return resonance( r, rho, g, penetrationShift, GX ); } );
    return lvalue::Type( std::move(rs), lstate.L(), offset );
  } else {
    auto rs =
      resonances( lstate.resonances(),
                  [&]( auto&& r )
                  { return resonance( r, rho, g, penetrationShift ); } );
    return lvalue::Type( std::move(rs), lstate.L() );
  }
}

template< typename ChannelRatio,
          typename StatisticalFactor >
static auto lvalue( const ENDF::resolved::SLBW::LState& lstate,
                    ChannelRatio&& rho,
                    StatisticalFactor&& g ){
  switch( lstate.L() ){
  case 0: return lvalue( lstate, rho, g, penetrationShift( Integer<0>{} ) );
  case 1: return lvalue( lstate, rho, g, penetrationShift( Integer<1>{} ) );
  case 2: return lvalue( lstate, rho, g, penetrationShift( Integer<2>{} ) );
  case 3: return lvalue( lstate, rho, g, penetrationShift( Integer<3>{} ) );
  case 4: return lvalue( lstate, rho, g, penetrationShift( Integer<4>{} ) );
  }
  throw std::exception();
}

template< typename BreitWigner,
          typename ChannelRatio,
          typename StatisticalFactor >
static auto lvalues( const BreitWigner& bw,
                     ChannelRatio&& rho,
                     StatisticalFactor&& g ){
  return
    bw.LStates()
    | ranges::view::transform
      ( [&]( auto&& lstate ){ return lvalue( lstate, rho, g ); } )
    | ranges::to_vector;
}
