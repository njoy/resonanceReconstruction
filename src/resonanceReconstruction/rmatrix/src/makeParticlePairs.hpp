unsigned int
findParticlePairForReaction(
    const endf::RMatrixLimited::ParticlePairs& endfPairs,
    int mt ) {

  const auto reactions = endfPairs.MT();
  auto found = std::find_if( ranges::cpp20::begin( reactions ),
                             ranges::cpp20::end( reactions ),
                             [&] ( const auto& reaction )
                                 { return mt == reaction; } );
  return std::distance( ranges::begin( reactions ), found );
}

unsigned int
incident( const endf::RMatrixLimited::ParticlePairs& endfPairs ) {

  // the incident particle pair is the one associated with mt2
  return findParticlePairForReaction( endfPairs, 2 );
}

unsigned int
eliminated( const endf::RMatrixLimited::ParticlePairs& endfPairs ) {

  // the incident particle pair is the one associated with mt102
  return findParticlePairForReaction( endfPairs, 102 );
}

std::vector< ParticlePair >
makeParticlePairs( const endf::RMatrixLimited::ParticlePairs& endfPairs,
                   const AtomicMass& neutronMass,
                   const ElectricalCharge& elementaryCharge,
                   const ParticleID& incident,
                   const ParticleID& target ) {

  // a few useful lambdas
  auto makeParticleIDs = [&] ( int mt )
    -> std::tuple< ParticleID, ParticleID, ParticlePairID > {

    elementary::ReactionType type( mt );
    if ( mt == 2 ) {

      return { incident, target, ParticlePairID( incident, target ) };
    }
    else if ( mt == 18 ) {

      return { incident, target, ParticlePairID( "fission" ) };
    }
    else {

      const auto particles = type.particles();
      if ( particles.size() > 1 ) {


      }
      auto outgoing = particles.front();
      auto residual = resolve( incident, target, type );
      return { outgoing, residual, ParticlePairID( outgoing, residual ) };
    }
  };
  auto first = [] ( const auto& tuple ) { return std::get< 0 >( tuple ); };
  auto second = [] ( const auto& tuple ) { return std::get< 1 >( tuple ); };
  auto third = [] ( const auto& tuple ) { return std::get< 2 >( tuple ); };
  auto makeParticle = [&] ( const ParticleID& id, double mass, double charge,
                            double spin, int parity ) {

    return Particle{ id, mass * neutronMass, charge * elementaryCharge,
                     std::abs( spin ),
                     spin == 0.0 ? ( parity >= 0 ? Parity( +1 ) : Parity( -1 ) )
                                 : ( spin > 0. ? Parity( +1 ) : Parity( -1 ) ) };
  };
  auto makeParticlePair = [&] ( const Particle& particle,
                                const Particle& residual,
                                const ParticlePairID& id ) {

    return ParticlePair{ particle, residual, id };
  };

  // do some range magic
  auto identifiers = endfPairs.MT() | ranges::views::transform( makeParticleIDs );
  auto particles = ranges::views::zip_with(
                       makeParticle,
                       identifiers | ranges::views::transform( first ),
                       endfPairs.massParticleA(),
                       endfPairs.chargeParticleA(),
                       endfPairs.spinParticleA(),
                       endfPairs.parityParticleA() );
  auto residuals = ranges::views::zip_with(
                       makeParticle,
                       identifiers | ranges::views::transform( second ),
                       endfPairs.massParticleB(),
                       endfPairs.chargeParticleB(),
                       endfPairs.spinParticleB(),
                       endfPairs.parityParticleB() );
  auto pairs = identifiers | ranges::views::transform( third );

  return ranges::to< std::vector< ParticlePair > >(
           ranges::views::zip_with(
             makeParticlePair,
             particles,
             residuals,
             pairs ) );
}
