template< typename Functor >
decltype(auto)
operator()( const ENDF::ResonanceRange& range,
            Functor&& callback ) const {

  try {
    const EnergyRange energyRange {
      range.EL() * electronVolts,
      range.EH() * electronVolts
    };
    const auto mlbw = std::get< ENDF::resolved::MLBW >( range.parameters() );

    /*
    Previously, the below code was structured like this:

      if( range.NRO() ){
        switch( range.NAPS() ){
          case 0: return callback( ... );
          case 1: return callback( ... );
          case 2: return callback( ... );
        }
      } else {
        switch( range.NAPS() ){
          case 0: return callback( ... );
          case 1: return callback( ... );
        }
      }

    We were able to modestly improve compile times by consolidating calls when
    we could. Within the present function, based at least on instantiations in
    the present test suite, these two expressions:

      channelRadius( mlbw.lValues().front().AWRI() )
      radius( mlbw.AP() )

    appear to always be of the same *type*. This allowed us to consolidate two
    calls, below, using the ternary conditional, with no change in meaning. We
    were unable to perform additional consolidations, as the types involved in
    doing so proved to be different, and thus incompatible with ?:.
    */

    if (range.NRO()) {
      switch (range.NAPS()) {
        case 0: // fall-through...
        case 2:
          return
            callback( build(
              energyRange, mlbw,
              range.NAPS() == 0
                ? channelRadius(mlbw.lValues().front().AWRI())
                : radius(mlbw.AP()),
              radius(range.scatteringRadius().value())
            ));
        case 1:
          return
            callback( build(
              energyRange, mlbw,
              radius(range.scatteringRadius().value())
            ));
      }
    } else {
      switch( range.NAPS() ){
        case 0:
          return
            callback( build(
              energyRange, mlbw,
              channelRadius(mlbw.lValues().front().AWRI()),
              radius(mlbw.AP())
            ));
        case 1:
          return
            callback( build(
              energyRange, mlbw,
              radius(mlbw.AP())
            ));
      }
    }

#if 0

    // ------------------------
    // Original code:
    // ------------------------

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

#endif

  } // try

  catch( std::bad_optional_access& ){
    Log::error( "Resonance range doesn't have scattering radius." );
    throw;
  }

  catch ( ... ) {
    throw std::runtime_error(
      "The resonance range does not appear to contain MLBW parameters" );
  }
}
