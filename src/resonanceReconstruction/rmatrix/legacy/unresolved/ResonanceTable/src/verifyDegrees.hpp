static
void verifyDegrees( const Degrees& degrees ) {

  // verify that the number of degrees is equal to the number of channels
  // verify that each channel is between 0 and 4

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
