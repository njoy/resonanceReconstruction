static
void verifyTable( const std::vector< Resonance >& resonances,
                  const Degrees& degrees ) {

  if ( resonances.size() < 2 ) {

    Log::error( "There must be at least two resonances in an unresolved "
                "resonance table to interpolate on the parameters" );
    Log::info( "Found {}", resonances.size() );
    throw std::exception();
  }

  const auto verifyDegreesOfFreedom = [] ( const auto& entry ) {

    if ( entry > 4 ) {

      Log::error( "The degrees of freedom for each channel must between 0 and 4" );
      Log::info( "Found {}", entry );
      throw std::exception();
    }
  };

  verifyDegreesOfFreedom( degrees.elastic );
  verifyDegreesOfFreedom( degrees.capture );
  verifyDegreesOfFreedom( degrees.fission );
  verifyDegreesOfFreedom( degrees.competition );
}
