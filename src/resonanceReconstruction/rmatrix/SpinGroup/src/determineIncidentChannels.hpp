static
std::vector< unsigned int >
determineIncidentChannels( const std::vector< ParticleChannel >& channels ) {

  std::vector< unsigned int > result;

  auto isIncidentChannel = [] ( const auto& channel )
                              { return channel.isIncidentChannel(); };

  auto iter = channels.begin();
  while ( ( iter =
              std::find_if(
                  iter, channels.end(),
                  [&] ( const ParticleChannel& channel )
                      { return std::visit( isIncidentChannel, channel ); } ) )
          != channels.end() ) {

    result.push_back( std::distance( channels.begin(), iter ) );
    ++iter;
  }
  return result;
}
