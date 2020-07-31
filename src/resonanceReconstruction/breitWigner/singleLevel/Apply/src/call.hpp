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

  } // try

  catch( std::bad_optional_access& ){
    Log::error( "Resonance range doesn't have scattering radius." );
    throw;
  }

  catch( ... ) {
    throw std::runtime_error(
      "The resonance range does not appear to contain SLBW parameters" );
  }
}
