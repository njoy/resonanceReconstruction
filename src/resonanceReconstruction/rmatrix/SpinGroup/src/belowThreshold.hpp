auto belowThreshold( const Energy& energy ) {

  Energy value = this->incidentPair().massRatio() * energy;
  return this->channels()
    | ranges::view::transform(
        [value] ( const auto& channel )
            { decltype(auto) q =
                std::visit( [] ( const auto& channel )
                               { return channel.particlePair().Q(); }, channel );
              return ( value + q ) < 0.0 * electronVolt; } );
}

