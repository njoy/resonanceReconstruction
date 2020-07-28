template< typename Functor >
decltype(auto)
operator()( const ENDF::ResonanceRange& range,
            Functor&& callback ) const {

  try {
    const EnergyRange energyRange {
      range.EL() * electronVolts,
      range.EH() * electronVolts
    };
    const auto slbw = std::get< ENDF::resolved::SLBW >( range.parameters() );

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

      channelRadius( slbw.lValues().front().AWRI() )
      radius( slbw.AP() )

    as well as these two:

      radius( range.scatteringRadius().value() )
      radius( slbw.AP() )

    appear to always be of the same *type*. This allowed us to consolidate the
    original five calls into two calls, below, using the ternary conditional,
    with no change in meaning.
    */

    return
      ( range.NRO() && range.NAPS() == 0) ||
      ( range.NRO() && range.NAPS() == 2) ||
      (!range.NRO() && range.NAPS() == 0)
    ? callback( build(
          energyRange, slbw,
          range.NAPS() == 0
            ? channelRadius( slbw.lValues().front().AWRI() )
            : radius( slbw.AP() ),
          range.NRO()
            ? radius( range.scatteringRadius().value() )
            : radius( slbw.AP() )
         )
       )
    : callback( build(
          energyRange, slbw,
          range.NRO()
            ? radius( range.scatteringRadius().value() )
            : radius( slbw.AP() )
          )
      );

#if 0

    // ------------------------
    // Original code:
    // ------------------------

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

#endif

  } // try

  catch( std::bad_optional_access& ){
    Log::error( "Resonance range doesn't have scattering radius." );
    throw;
  }

  catch ( ... ) {
    throw std::runtime_error(
      "The resonance range does not appear to contain SLBW parameters" );
  }
}
