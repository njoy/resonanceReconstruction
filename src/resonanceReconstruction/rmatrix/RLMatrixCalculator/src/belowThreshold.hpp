template < typename Channels >
auto belowThreshold( const Energy& energy, Channels channels ) {

  auto incident = [] ( const auto& channel )
                     { return channel.incident(); };
  auto particlePair = [] ( const auto& channel )
                         { return channel.particlePair(); };

  unsigned int index =
       ranges::distance(
           ranges::begin( channels ),
           ranges::find_if(
               channels,
               [&] ( const auto& channel )
                   { return std::visit( incident, channel ); } ) );

  auto incidentPair = std::visit( particlePair, channels[index] );
  Energy value = incidentPair.massRatio() * energy;
  return channels
    | ranges::view::transform(
        [value] ( const auto& channel )
            { decltype(auto) q =
                std::visit( [] ( const auto& channel )
                               { return channel.Q(); }, channel );
              return ( value + q ) < 0.0 * electronVolt; } );
}
