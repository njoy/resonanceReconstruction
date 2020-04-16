template< typename Functor >
decltype(auto)
operator()( const ENDF::ResonanceRange& range,
            Functor&& callback ) const {
  try {

    EnergyRange energyRange{ range.EL() * electronVolts,
                            range.EH() * electronVolts };
    auto mlbw = std::get< ENDF::resolved::MLBW >( range.parameters() );

    if( range.NRO() ){
      switch( range.NAPS() ){
      case 0:
        return callback( build( energyRange, mlbw,
                                channelRadius( mlbw.lValues().front().AWRI() ),
                                radius( mlbw.AP() ) ) );
      case 1:
        return callback( build( energyRange, mlbw, radius( mlbw.AP() ) ) );
      }
    } else {
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
    }
  }
  catch ( ... ) {

    Log::error( 
      "The resonance range does not appear to contain MLBW parameters" );
    throw;
  }
}
