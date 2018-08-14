static
void verifyQuantumNumbers( const std::vector< ParticleChannel >& channels ) {

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

