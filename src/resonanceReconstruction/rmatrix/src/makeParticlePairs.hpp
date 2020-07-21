inline unsigned int
findParticlePairForReaction(
    const ENDF::resolved::RMatrixLimited::ParticlePairs& endfPairs,
    int mt ) {

  const auto reactions = endfPairs.MT();
  auto found = std::find_if( ranges::begin( reactions ),
                             ranges::end( reactions ),
                             [&] ( const auto& reaction )
                                 { return mt == reaction; } );
  return std::distance( ranges::begin( reactions ), found );
}

inline unsigned int
incident( const ENDF::resolved::RMatrixLimited::ParticlePairs& endfPairs ) {

  // the incident particle pair is the one associated with mt2
  return findParticlePairForReaction( endfPairs, 2 );
}

inline unsigned int
eliminated( const ENDF::resolved::RMatrixLimited::ParticlePairs& endfPairs ) {

  // the incident particle pair is the one associated with mt102
  return findParticlePairForReaction( endfPairs, 102 );
}

inline std::vector< ParticlePair >
makeParticlePairs( const ENDF::resolved::RMatrixLimited::ParticlePairs& endfPairs,
                   const AtomicMass& neutronMass,
                   const ElectricalCharge& elementaryCharge ) {

  const std::map< int, std::string > particleID = {
    { 0, "g" }, { 1, "n" }, { 1001, "p" }, { 1002, "h2" }, { 1003, "h3" },
    { 2003, "he3" }, { 2004, "he4" } };

  // a few useful lambdas
  auto makeParticleIDs = [&] ( double za, double ma,
                               double zb, double mb,
                               int mt )  -> std::pair< ParticleID, ParticleID > {

    unsigned int aa = std::round( AtomicMass( ma * neutronMass ).value );
    unsigned int ab = std::round( AtomicMass( mb * neutronMass ).value );

    std::string level = "_e0";
    if ( mt == 51 )  { level = "_e1"; }
    else if ( mt == 600 ) { level = "_e1"; };

    return { ParticleID(
               particleID.at( static_cast< unsigned int >( za ) * 1000 + aa ) ),
             ParticleID(
               elementary::ElementID( static_cast< unsigned int >( zb ) ).symbol()
               + std::to_string( ab )
               + level ) };
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
                                const Particle& residual ) {

    return ParticlePair{ particle, residual };
  };

  // do some range magic
  auto identifiers = ranges::view::zip_with(
                         makeParticleIDs,
                         endfPairs.chargeParticleA(),
                         endfPairs.massParticleA(),
                         endfPairs.chargeParticleB(),
                         endfPairs.massParticleB(),
                         endfPairs.MT() );
  auto particles = ranges::view::zip_with(
                       makeParticle,
                       identifiers | ranges::view::transform( first ),
                       endfPairs.massParticleA(),
                       endfPairs.chargeParticleA(),
                       endfPairs.spinParticleA(),
                       endfPairs.parityParticleA() );
  auto residuals = ranges::view::zip_with(
                       makeParticle,
                       identifiers | ranges::view::transform( second ),
                       endfPairs.massParticleB(),
                       endfPairs.chargeParticleB(),
                       endfPairs.spinParticleB(),
                       endfPairs.parityParticleB() );

  return ranges::view::zip_with(
             makeParticlePair,
             particles,
             residuals );
}
