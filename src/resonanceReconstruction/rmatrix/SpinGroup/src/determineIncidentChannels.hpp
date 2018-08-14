static
std::vector< unsigned int >
determineIncidentChannels( const ParticlePairID& incident,
                           const std::vector< ParticleChannel >& channels ) {

  std::vector< unsigned int > result;

  auto iter = channels.begin();
  while ( ( iter =
              std::find_if(
                  iter, channels.end(),
                  [&] ( const ParticleChannel& input )
                      { return
                            std::visit(
                                [&] ( const auto& channel )
                                    { return channel.particlePair().pairID() ==
                                             incident; }, input ); } ) ) != 
          channels.end() )
  {
    result.push_back( std::distance( channels.begin(), iter ) );
    ++iter;
  }
  return result;
}

