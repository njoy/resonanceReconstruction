template< typename Functor >
decltype(auto)
operator()( const ENDF::ResonanceRange& range,
            Functor&& callback ) const {

  try {

    EnergyRange energyRange{ range.EL() * electronVolts,
                            range.EH() * electronVolts };
    auto slbw = std::get< ENDF::resolved::SLBW >( range.parameters() );

    switch( range.NRO() ){
    case 0:
      switch( range.NAPS() ){
      case 0:
        return callback( build( energyRange, slbw,
                                channelRadius( slbw.lValues().front().AWRI() ),
                                radius( slbw.AP() ) ) );
      case 1:
        return callback( build( energyRange, slbw, radius( slbw.AP() ) ) );
      }
    case 1:
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
    }
  }
  catch ( ... ) {

    Log::error( "The resonance range does not appear to contain SLBW parameters" );
    throw std::exception();
  }
}
