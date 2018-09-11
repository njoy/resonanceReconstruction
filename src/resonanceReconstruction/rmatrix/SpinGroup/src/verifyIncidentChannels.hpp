static
void verifyIncidentChannels( const std::vector< ParticleChannel >& channels,
                             const std::vector< unsigned int >& indices ) {

  // verify that there is at least one incident channel and not more than
  // the actual number of channels

  const unsigned int size = channels.size();
  const unsigned int number = indices.size();
  if ( number == 0 ) {

    Log::error( "No incident channels are given or could be found." );
    throw std::exception();
  }
  if ( number > size ) {

    Log::error( "More incident channels than there are actual channels." );
    Log::info( "Expected {} or less indices, found {} indices", size, number );
    throw std::exception();
  }

  // verify the indices

  auto checkIndex = [&] ( const unsigned int index ) {

    if ( index > size - 1 ) {

      Log::error( "Erroneous incident channel index." );
      Log::info( "Expected an index less than or equal to {}, found {}",
                 size - 1, index );
      throw std::exception();
    }
  };

  ranges::for_each( indices, checkIndex );

  // verify that each incident channel has the same incident pair

  const auto getParticlePairID = [] ( const auto& entry ) {

    return std::visit( [] ( const auto& channel )
                          { return channel.particlePair().pairID(); },
                       entry );
  };

  const auto reference = getParticlePairID( channels[ indices.front() ] );

  auto checkIncidentPairID = [&] ( const auto& entry ) {

    const auto current = getParticlePairID( entry );
    if ( current != reference ) {

      Log::error( "Inconsistent incident particle pair." );
      Log::info( "Expected particle pair {}, found {}", reference, current );
      throw std::exception();
    }
  };

  ranges::for_each( indices | ranges::view::drop_exactly( 1 )
                            | ranges::view::transform(
                                  [&] ( const unsigned int index )
                                      { return channels[ index ]; } ),
                    checkIncidentPairID );
}
