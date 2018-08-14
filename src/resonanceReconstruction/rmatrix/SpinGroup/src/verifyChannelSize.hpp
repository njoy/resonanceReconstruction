static
void verifyChannelSize( const unsigned int size ) {

  if ( size == 0 ) {

    Log::error( "The number of channels in a spin group cannot be 0." );
    throw std::exception();
  }
}

