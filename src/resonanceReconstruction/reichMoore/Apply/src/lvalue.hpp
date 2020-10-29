template< typename WaveNumber,
          typename Radius,
          typename StatisticalFactor,
          typename PenetrationShift >
static auto lvalue( const endf::ReichMoore::LValue& lstate,
                    WaveNumber&& k,
                    Radius&& r,
                    StatisticalFactor&& g,
                    PenetrationShift&& penetrationShift,
                    bool useAPL ){
    auto rs =
      ( lstate.APL() != 0.0 and useAPL ) ?
      [&]{
      auto rho = [ &, lr = radius( lstate.APL() ) ]
                   ( auto&& energy ){ return lr(energy) * k( energy ); };
        return resonances
               ( lstate.resonances(),
                 [ & ]( auto&& res )
                 { return resonance( res, rho, g, penetrationShift ); } );
      }() :
      [&]{
        auto rho = [&]( auto&& energy ){ return r(energy) * k( energy ); };
        return resonances
               ( lstate.resonances(),
                 [ & ]( auto&& res )
                 { return resonance( res, rho, g, penetrationShift ); } );
      }();

    return lstate.APL() != 0.0 ?
      Lvalue( std::move(rs), lstate.L(), radius( lstate.APL() ) ) :
      Lvalue( std::move(rs), lstate.L() );
}

template< typename WaveNumber,
          typename Radius,
          typename StatisticalFactor >
static auto lvalue( const endf::ReichMoore::LValue& lstate,
                    WaveNumber&& k,
                    Radius&& r,
                    StatisticalFactor&& g,
                    bool useAPL ){
  switch( lstate.L() ){
  case 0:
    return lvalue( lstate, k, r, g, penetrationShift( Integer<0>{} ), useAPL );
  case 1:
    return lvalue( lstate, k, r, g, penetrationShift( Integer<1>{} ), useAPL );
  case 2:
    return lvalue( lstate, k, r, g, penetrationShift( Integer<2>{} ), useAPL );
  case 3:
    return lvalue( lstate, k, r, g, penetrationShift( Integer<3>{} ), useAPL );
  case 4:
    return lvalue( lstate, k, r, g, penetrationShift( Integer<4>{} ), useAPL );
  }
  throw std::exception();
}

template< typename WaveNumber,
          typename Radius,
          typename StatisticalFactor >
static auto lvalues( const endf::ReichMoore& rm,
                     WaveNumber&& k,
                     Radius&& r,
                     StatisticalFactor&& g,
                     bool useAPL ){
  return
    rm.lValues()
    | ranges::view::transform
    ( [&]( auto&& lstate ){ return lvalue( lstate, k, r, g, useAPL ); } )
    | ranges::to_vector;
}
