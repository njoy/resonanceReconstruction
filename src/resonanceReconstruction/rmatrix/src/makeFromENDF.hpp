inline unsigned int
incident( const ENDF::resolved::RMatrixLimited::ParticlePairs& pairs ) {

  // determine the incident particle pair (the particle pair with mt == 2)
  const auto reactions = pairs.MT();
  auto found = std::find_if( ranges::begin( reactions ),
                             ranges::end( reactions ),
                             [&] ( const auto& mt ) { return mt == 2; } );
  return std::distance( ranges::begin( reactions ), found );
}

inline std::vector< ParticlePair >
makeParticlePairs( const ENDF::resolved::RMatrixLimited::ParticlePairs& pairs,
                   const AtomicMass& neutronMass,
                   const ElectricalCharge& elementaryCharge ) {

  // a few useful lambdas
  auto makeParticleIdentifiers = [] ( const auto& pairs ) {

    //! @todo replace these temporary identifiers
    return pairs.MT()
             | ranges::view::transform(
                   [] ( int mt ) -> std::pair< ParticleID, ParticleID > {

                     std::string id = "mt" + std::to_string( mt );
                     return { id + "-a", id + "-b" };
                   } );
  };
  auto first = [] ( const auto& pair ) { return pair.first; };
  auto second = [] ( const auto& pair ) { return pair.second; };
  auto makeParticle = [&] ( const ParticleID& id, double mass, double charge,
                            double spin, int parity ) {

    return Particle{ id, mass * neutronMass, charge * elementaryCharge,
                     std::abs( spin ),
                     spin == 0.0 ? ( parity >= 0 ? Parity( +1 ) : Parity( -1 ) )
                                 : ( spin > 0. ? Parity( +1 ) : Parity( -1 ) ) };
  };
  auto makeParticlePair = [&] ( const Particle& particle,
                                const Particle& residual,
                                double Q ) {

    return ParticlePair{ particle, residual, Q * electronVolt };
  };

  // do some range magic
  auto identifiers = makeParticleIdentifiers( pairs );
  auto particles = ranges::view::zip_with(
                       makeParticle,
                       identifiers | ranges::view::transform( first ),
                       pairs.massParticleA(),
                       pairs.chargeParticleA(),
                       pairs.spinParticleA(),
                       pairs.parityParticleA() );
  auto residuals = ranges::view::zip_with(
                       makeParticle,
                       identifiers | ranges::view::transform( second ),
                       pairs.massParticleB(),
                       pairs.chargeParticleB(),
                       pairs.spinParticleB(),
                       pairs.parityParticleB() );

  return ranges::view::zip_with(
             makeParticlePair,
             particles,
             residuals,
             pairs.Q() );
}
