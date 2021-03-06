template< typename Functor >
decltype(auto)
operator()( const endf::ResonanceRange& range,
            Functor&& callback ) const {

  try {

    EnergyRange energyRange{ range.EL() * electronVolts,
                            range.EH() * electronVolts };
    auto slbw = std::get< endf::SingleLevelBreitWigner >( range.parameters() );

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
