static
std::vector< unsigned int >
determineIncidentChannels( const ParticlePairID& incident,
                           std::vector< ParticleChannel >& channels ) {

  std::vector< unsigned int > result;

  for ( ParticleChannel& channel : channels ) {

    std::visit( [] ( auto& channel )
                   { if ( channel.incident() ) { channel.toggleIncident(); } },
                channel );
  }

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
          channels.end() ) {

    result.push_back( std::distance( channels.begin(), iter ) );
    std::visit( [] ( auto& channel ) { channel.toggleIncident(); }, *iter );
    ++iter;
  }
  return result;
}
