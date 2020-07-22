template< typename Functor >
decltype(auto)
operator()( const ENDF::ResonanceRange& range,
            Functor&& callback ) const {

  try {

    EnergyRange energyRange{ range.EL() * electronVolts,
                             range.EH() * electronVolts };
    auto rm = std::get< ENDF::resolved::RM >( range.parameters() );

#if 0

    /*
    // This works:
    decltype(auto) foo = condition
      ? channelRadius( rm.lValues().front().AWRI() )
      : radius( rm.AP() );
    // This FAILS:
    decltype(auto) bar = condition
      ? radius( range.scatteringRadius().value() )
      : radius( rm.AP() );
    */

    decltype(auto) foo = channelRadius( rm.lValues().front().AWRI() );
    decltype(auto) bar = radius( range.scatteringRadius().value() );
    decltype(auto) baz = radius( rm.AP() );

    if( range.NRO() )
      switch( range.NAPS() ){
      case 0: return callback(build(energyRange,rm,Neither{},foo,bar,false));
      case 1: return callback(build(energyRange,rm,Neither{},bar,    false));
      case 2: return callback(build(energyRange,rm,Channel{},baz,bar,true ));
      }
    else
      switch( range.NAPS() ){
      case 0: return callback(build(energyRange,rm,Scattering{},foo,baz,false));
      case 1: return callback(build(energyRange,rm,Both{},baz,true));
      }

#endif

    if( range.NRO() ){
      switch( range.NAPS() ){
        case 0:
          return callback( build( energyRange,
                                  rm,
                                  Neither{},
                                  channelRadius( rm.lValues().front().AWRI() ),
                                  radius( range.scatteringRadius().value() ),
                                  false ) );
        case 1:
          return callback( build( energyRange,
                                  rm,
                                  Neither{},
                                  radius( range.scatteringRadius().value() ),
                                  false ) );
        case 2:
          return callback( build( energyRange,
                                  rm,
                                  Channel{},
                                  radius( rm.AP() ),
                                  radius( range.scatteringRadius().value() ),
                                  true ) );
      }
    } else {
      switch( range.NAPS() ){
        case 0:
          return callback( build( energyRange,
                                  rm,
                                  Scattering{},
                                  channelRadius( rm.lValues().front().AWRI() ),
                                  radius( rm.AP() ),
                                  false ) );
        case 1:
          return callback( build( energyRange,
                                  rm, Both{}, radius( rm.AP() ), true ) );
      }
    }
  }

  catch ( ... ) {

    throw std::runtime_error( 
      "The resonance range does not appear to contain Reich-Moore parameters" );
  }
}
