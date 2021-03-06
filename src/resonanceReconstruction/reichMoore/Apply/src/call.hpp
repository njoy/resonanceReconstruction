template< typename Functor >
decltype(auto)
operator()( const endf::ResonanceRange& range,
            Functor&& callback ) const {

  try {

    EnergyRange energyRange{ range.EL() * electronVolts,
                             range.EH() * electronVolts };
    auto rm = std::get< endf::ReichMoore >( range.parameters() );

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
