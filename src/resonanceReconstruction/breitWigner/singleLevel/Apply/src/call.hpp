template< typename Functor >
decltype(auto)
operator()( const ENDF::ResonanceRange& range,
            Functor&& callback ) const {

  try {

    EnergyRange energyRange{ range.EL() * electronVolts,
                            range.EH() * electronVolts };
    auto slbw = std::get< ENDF::resolved::SLBW >( range.parameters() );

    /*
    // These work:
    decltype(auto) foo = condition
      ? channelRadius( slbw.lValues().front().AWRI() )
      : radius( slbw.AP() );
    decltype(auto) bar = condition
      ? radius( range.scatteringRadius().value() )
      : radius( slbw.AP() );
    */

#define foo channelRadius( slbw.lValues().front().AWRI() )
#define bar radius( range.scatteringRadius().value() )
#define baz radius( slbw.AP() )

    if (( range.NRO() && range.NAPS() == 0) ||
        ( range.NRO() && range.NAPS() == 2) ||
        (!range.NRO() && range.NAPS() == 0))
       return callback( build( energyRange, slbw, range.NAPS() == 0 ? foo : baz, range.NRO() ? bar : baz ) );
    else
       return callback( build( energyRange, slbw, range.NRO() ? bar : baz ) );

#undef foo
#undef bar
#undef baz

    /*
    if ( range.NRO() && range.NAPS() == 0)
       return callback( build( energyRange, slbw, foo, bar ) );
    if ( range.NRO() && range.NAPS() == 2)
       return callback( build( energyRange, slbw, baz, bar ) );
    if (!range.NRO() && range.NAPS() == 0)
       return callback( build( energyRange, slbw, foo, baz ) );
    if ( range.NRO() && range.NAPS() == 1)
       return callback( build( energyRange, slbw, bar ) );
    if (!range.NRO() && range.NAPS() == 1)
       return callback( build( energyRange, slbw, baz ) );
    */

    /*
    if( range.NRO() )
      switch( range.NAPS() ){
        case 0: return callback( build( energyRange, slbw, foo, bar ) );
        case 1: return callback( build( energyRange, slbw, bar ) );
        case 2: return callback( build( energyRange, slbw, baz, bar ) );
      }
    else
      switch( range.NAPS() ){
        case 0: return callback( build( energyRange, slbw, foo, baz ) );
        case 1: return callback( build( energyRange, slbw, baz ) );
      }
    */

#if 0
    decltype(auto) foo = channelRadius( slbw.lValues().front().AWRI() );
    decltype(auto) bar = radius( range.scatteringRadius().value() );
    decltype(auto) baz = radius( slbw.AP() );

    if( range.NRO() )
      switch( range.NAPS() ){
        case 0: return callback( build( energyRange, slbw, foo, bar ) );
        case 1: return callback( build( energyRange, slbw, bar ) );
        case 2: return callback( build( energyRange, slbw, baz, bar ) );
      }
    else
      switch( range.NAPS() ){
        case 0: return callback( build( energyRange, slbw, foo, baz ) );
        case 1: return callback( build( energyRange, slbw, baz ) );
      }
#endif

    /*
    if( range.NRO() ){
      switch( range.NAPS() ){
      case 0:
        return callback( build( energyRange, slbw,
                                channelRadius( slbw.lValues().front().AWRI() ),
                                radius( range.scatteringRadius().value() ) ) );
      case 1:
        return callback( build( energyRange, slbw,
                                radius( range.scatteringRadius().value() ) ) );
      case 2:
        return callback( build( energyRange, slbw, radius( slbw.AP() ),
                                radius( range.scatteringRadius().value() ) ) );
      }
    } else {
      switch( range.NAPS() ){
      case 0:
        return callback( build( energyRange, slbw,
                                channelRadius( slbw.lValues().front().AWRI() ),
                                radius( slbw.AP() ) ) );
      case 1:
        return callback( build( energyRange, slbw, radius( slbw.AP() ) ) );
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
      "The resonance range does not appear to contain SLBW parameters" );
  }
}
