static
std::vector< unsigned int >
determineIncidentChannels( const std::vector< ParticleChannel >& channels ) {

  std::vector< unsigned int > result;

  const auto isIncidentChannel = [] ( const auto& channel ) {

    return channel.isIncidentChannel();
  };

  const auto getIsIncidentChannel = [&] ( const auto& channel ) {

    return std::visit( isIncidentChannel, channel );
  };

  auto iter = channels.begin();
  while ( ( iter = std::find_if( iter, channels.end(), getIsIncidentChannel ) )
          != channels.end() ) {

    result.push_back( std::distance( channels.begin(), iter ) );
    ++iter;
  }
  return result;
}
