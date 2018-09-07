static
void verifyChannels( const std::vector< ParticleChannel >& channels ) {

  // verify that there is at least one channel

  if ( channels.size() == 0 ) {

    Log::error( "The number of channels in a spin group cannot be 0." );
    throw std::exception();
  }

  // verify that each channel is unique

  const auto verifyUniqueChannel = [&] ( const auto& entry ) {
      
    const auto getChannelID = [] ( const auto& entry ) {

      return std::visit( [] ( const auto& channel )
                            { return channel.channelID(); },
                         entry );
    };

    const auto label = getChannelID( entry );
    if ( ranges::count_if(
             channels,
             [&] ( const auto& channel )
                 { return getChannelID( channel ) == label; } ) > 1 ) {

      Log::error( "Channels in the spin group do not seem to be unique." );
      Log::info( "Channel {} is present at least twice", label );
      throw std::exception();
    }
  };

  ranges::for_each( channels, verifyUniqueChannel );

  // verify that each channel's Jpi is the same

  const auto getQuantumNumbers = [] ( const auto& entry ) {
    return std::visit( [] ( const auto& channel )
                          { return channel.quantumNumbers(); },
                       entry );
  };

  const auto reference = getQuantumNumbers( channels.front() );

  auto checkQuantumNumbers = [&] ( const auto& entry ) {
    const auto current = getQuantumNumbers( entry );
    const bool mismatch = 
      ( ( current.totalAngularMomentum() != reference.totalAngularMomentum() ) ||
        ( current.parity() != reference.parity() ) );

    if ( mismatch ) {

      Log::error( "Total angular momentum and parity mismatch in spin group." );
      Log::info( "Total angular momentum J: expected {}, found {}",
                 reference.totalAngularMomentum(),
                 current.totalAngularMomentum() );
      Log::info( "Parity: expected {}, found {}",
                 reference.parity(),
                 current.parity() );
      throw std::exception();
    }
  };

  ranges::for_each( channels | ranges::view::drop_exactly( 1 ),
                    checkQuantumNumbers );
}
