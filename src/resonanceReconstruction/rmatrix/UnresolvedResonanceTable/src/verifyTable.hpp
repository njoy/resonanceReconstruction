static
void verifyTable( const std::vector< ChannelID >& channels,
                  const std::vector< unsigned int >& degrees ) {

  // verify that the number of degrees is equal to the number of channels
  // verify that each channel is between 0 and 4

  unsigned int expected = channels.size();
  unsigned int number = degrees.size();
  if ( expected != number ) {

    Log::error( "Number of degrees of freedom different from expected." );
    Log::info( "Found {}, expected {}", number, expected );
    throw std::exception();
  }

  const auto verifyDegreesOfFreedom = [] ( const auto& entry ) {

    if ( entry > 4 ) {

      Log::error( "The degrees of freedom for each channel must between 0 and 4" );
      Log::info( "Found {}", entry );
      throw std::exception();
    }
  };

  ranges::for_each( degrees, verifyDegreesOfFreedom );
}
