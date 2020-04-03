template< typename Functor >
decltype(auto)
operator()( const ENDF::ResonanceRange& range,
            Functor&& callback ) const {
  try {

    EnergyRange energyRange{ range.EL() * electronVolts,
                            range.EH() * electronVolts };
    auto mlbw = std::get< ENDF::resolved::MLBW >( range.parameters() );

    switch( range.NRO() ){
    case 0:
      switch( range.NAPS() ){
      case 0:
        return callback( build( energyRange, mlbw,
                                channelRadius( mlbw.lValues().front().AWRI() ),
                                radius( mlbw.AP() ) ) );
      case 1:
        return callback( build( energyRange, mlbw, radius( mlbw.AP() ) ) );
      }
    case 1:
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

    throw std::runtime_error( "The resonance range does not appear to contain "
                              "MLBW parameters" );
  }
}
