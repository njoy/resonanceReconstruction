static void
verifyIncidentParticlePairs( const std::vector< ParticleChannel >& channels ) {

  const auto incidentPairID = [] ( const auto& channel ) -> decltype(auto) {

    return channel.incidentParticlePair().pairID();
  };

  const auto getIncidentPairID = [&] ( const auto& channel ) -> decltype(auto) {

    return std::visit( incidentPairID, channel );
  };

  decltype(auto) in = getIncidentPairID( channels.front() );

  const auto equalIncidentPairID = [&] ( const auto& channel ) {

    return getIncidentPairID( channel ) == in;
  };

  unsigned int number = ranges::cpp20::count_if( channels, equalIncidentPairID );
  if ( number != channels.size() ) {

    Log::error( "Not all channels in the spin group have the same incident "
                "particle pair." );
    throw std::exception();
  }
}

static void
verifyUniqueChannels( const std::vector< ParticleChannel >& channels ) {

  const auto verifyUniqueChannel = [&] ( const auto& entry ) {

    const auto channelID = [] ( const auto& channel ) -> decltype(auto) {

      return channel.channelID();
    };

    const auto getChannelID = [&] ( const auto& channel ) -> decltype(auto) {

      return std::visit( channelID, channel );
    };

    const auto label = getChannelID( entry );

    const auto equalChannelID = [&] ( const auto& channel ) {

      return getChannelID( channel ) == label;
    };

    if ( ranges::cpp20::count_if( channels, equalChannelID ) > 1 ) {

      Log::error( "Channels in the spin group do not seem to be unique." );
      Log::info( "Channel {} is present at least twice", label );
      throw std::exception();
    }
  };

  ranges::cpp20::for_each( channels, verifyUniqueChannel );
}

static void
verifySpinParity( const std::vector< ParticleChannel >& channels ) {

  const auto quantumNumbers = [] ( const auto& channel ) -> decltype(auto) {

    return channel.quantumNumbers();
  };

  const auto getQuantumNumbers = [&] ( const auto& channel ) -> decltype(auto) {

    return std::visit( quantumNumbers, channel );
  };

  decltype(auto) reference = getQuantumNumbers( channels.front() );

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

  ranges::cpp20::for_each( channels, checkQuantumNumbers );
}

static
void verifyChannels( const std::vector< ParticleChannel >& channels ) {

  // verify that there is at least one channel
  if ( channels.size() == 0 ) {

    Log::error( "The number of channels in a spin group cannot be 0." );
    throw std::exception();
  }

  // verify that each channel has the same incident particle pair
  verifyIncidentParticlePairs( channels );

  // verify that each channel is unique
  verifyUniqueChannels( channels );

  // verify that each channel's Jpi is the same
  verifySpinParity( channels );
}
