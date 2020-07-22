template< typename Functor >
decltype(auto)
operator()( const ENDF::ResonanceRange& range,
            Functor&& callback ) const {

  try {

    EnergyRange energyRange{ range.EL() * electronVolts,
                            range.EH() * electronVolts };
    auto mlbw = std::get< ENDF::resolved::MLBW >( range.parameters() );

    if (range.NRO())
      switch (range.NAPS()) {
        case 0:
        case 2:
          return
            callback(
              build(
                energyRange,
                mlbw,
                range.NAPS() == 0
                  ? channelRadius(mlbw.lValues().front().AWRI())
                  : radius(mlbw.AP()),
                radius(range.scatteringRadius().value())
              )
            );

        case 1:
          return
            callback(
              build(
                energyRange,
                mlbw,
                radius(range.scatteringRadius().value())
              )
            );
      }
    else
      switch( range.NAPS() ){
        case 0:
          return
            callback(
              build(
                energyRange,
                mlbw,
                channelRadius(mlbw.lValues().front().AWRI()),
                radius(mlbw.AP())
              )
            );
        case 1:
          return
            callback(
              build(
                energyRange,
                mlbw,
                radius(mlbw.AP())
              )
            );
      }

#if 0

    /*
    // This works:
    decltype(auto) foo = condition
      ? channelRadius( mlbw.lValues().front().AWRI() )
      : radius( mlbw.AP() );
    // This FAILS:
    decltype(auto) bar = condition
      ? radius( range.scatteringRadius().value() )
      : radius( mlbw.AP() );
    */

    decltype(auto) foo = channelRadius( mlbw.lValues().front().AWRI() );
    decltype(auto) bar = radius( range.scatteringRadius().value() );
    decltype(auto) baz = radius( mlbw.AP() );

    if( range.NRO() )
      switch( range.NAPS() ){
        case 0: return callback( build( energyRange, mlbw, foo, bar ) );
        case 1: return callback( build( energyRange, mlbw, bar ) );
        case 2: return callback( build( energyRange, mlbw, baz, bar ) );
      }
    else
      switch( range.NAPS() ){
        case 0: return callback( build( energyRange, mlbw, foo, baz ) );
        case 1: return callback( build( energyRange, mlbw, baz ) );
      }

#endif

    /*
    if( range.NRO() ){
      switch( range.NAPS() ){
      case 0:
        return callback( build( energyRange, mlbw,
                                channelRadius( mlbw.lValues().front().AWRI() ),
                                radius( range.scatteringRadius().value() ) ) );
      case 1:
        return callback( build( energyRange, mlbw,
                                radius( range.scatteringRadius().value() ) ) );
      case 2:
        return callback( build( energyRange, mlbw, radius( mlbw.AP() ),
                                radius( range.scatteringRadius().value() ) ) );
      }
    } else {
      switch( range.NAPS() ){
      case 0:
        return callback( build( energyRange, mlbw,
                                channelRadius( mlbw.lValues().front().AWRI() ),
                                radius( mlbw.AP() ) ) );
      case 1:
        return callback( build( energyRange, mlbw, radius( mlbw.AP() ) ) );
      }
    }
    */
  }

  catch( std::bad_optional_access& ){
    Log::error( "Resonance range doesn't have scattering radius." );
    throw;
  }

  catch ( ... ) {
    throw std::runtime_error( 
      "The resonance range does not appear to contain MLBW parameters" );
  }
}
