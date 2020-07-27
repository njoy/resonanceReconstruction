static
void verifyIncidentChannels( const std::vector< unsigned int >& indices ) {

  // verify that there is at least one incident channel
  // indices are determined automatically, no need to test anything else

  const unsigned int number = indices.size();
  if ( number == 0 ) {

    Log::error( "No incident channels are given or could be found." );
    throw std::exception();
  }
}
