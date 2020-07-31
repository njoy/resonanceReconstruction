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

  } // try

  catch( std::bad_optional_access& ){
    Log::error( "Resonance range doesn't have scattering radius." );
    throw;
  }

  catch( ... ) {
    throw std::runtime_error(
      "The resonance range does not appear to contain MLBW parameters" );
  }
}
