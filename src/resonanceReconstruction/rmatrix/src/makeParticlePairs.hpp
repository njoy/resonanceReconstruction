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
    { 0, "g" }, { 1, "n" }, { 1001, "p" }, { 1002, "d" }, { 1003, "t" },
    { 2003, "he3" }, { 2004, "a" } };
  const std::map< int, std::string > elementID = {
    { 1, "H" }, { 2, "He" }, { 3, "Li" }, { 4, "Be" }, { 5, "B" },
    { 6, "C" }, { 7, "N" }, { 8, "O" }, { 9, "F" }, { 10, "Ne" },
    { 11, "Na" }, { 12, "Mg" }, { 13, "Al" }, { 14, "Si" }, { 15, "P" },
    { 16, "S" }, { 17, "Cl" }, { 18, "Ar" }, { 19, "K" }, { 20, "Ca" },
    { 21, "Sc" }, { 22, "Ti" }, { 23, "V" }, { 24, "Cr" }, { 25, "Mn" },
    { 26, "Fe" }, { 27, "Co" }, { 28, "Ni" }, { 29, "Cu" }, { 30, "Zn" },
    { 31, "Ga" }, { 32, "Ge" }, { 33, "As" }, { 34, "Se" }, { 35, "Br" },
    { 36, "Kr" }, { 37, "Rb" }, { 38, "Sr" }, { 39, "Y" }, { 41, "Nb" },
    { 40, "Zr" }, { 42, "Mo" }, { 43, "Tc" }, { 44, "Ru" }, { 45, "Rh" },
    { 46, "Pd" }, { 47, "Ag" }, { 48, "Cd" }, { 49, "In" }, { 50, "Sn" },
    { 51, "Sb" }, { 52, "Te" }, { 53, "I" }, { 55, "Cs" }, { 54, "Xe" },
    { 56, "Ba" }, { 57, "La" }, { 58, "Ce" }, { 60, "Nd" }, { 59, "Pr" },
    { 61, "Pm" }, { 62, "Sm" }, { 63, "Eu" }, { 64, "Gd" }, { 65, "Tb" },
    { 66, "Dy" }, { 67, "Ho" }, { 68, "Er" }, { 69, "Tm" }, { 70, "Yb" },
    { 71, "Lu" }, { 72, "Hf" }, { 73, "Ta" }, { 74, "W" }, { 75, "Re" },
    { 76, "Os" }, { 77, "Ir" }, { 78, "Pt" }, { 79, "Au" }, { 81, "Tl" },
    { 80, "Hg" }, { 82, "Pb" }, { 83, "Bi" }, { 84, "Po" }, { 85, "At" },
    { 86, "Rn" }, { 87, "Fr" }, { 88, "Ra" }, { 89, "Ac" }, { 91, "Pa" },
    { 90, "Th" }, { 92, "U" }, { 93, "Np" }, { 94, "Pu" }, { 95, "Am" },
    { 96, "Cm" }, { 97, "Bk" }, { 98, "Cf" }, { 99, "Es" }, { 100, "Fm" },
    { 101, "Md" }, { 102, "No" }, { 104, "Rf" }, { 103, "Lr" }, { 105, "Db" },
    { 106, "Sg" }, { 108, "Hs" }, { 107, "Bh" }, { 109, "Mt" }, { 110, "Ds" },
    { 111, "Rg" }, { 113, "Nh" }, { 112, "Cn" }, { 115, "Mc" }, { 114, "Fl" },
    { 116, "Lv" }, { 117, "Ts" }, { 118, "Og" }, { 119, "Uue" }, { 120, "Ubn" } };

  // a few useful lambdas
  auto makeParticleIDs = [&] ( double za, double ma,
                               double zb, double mb,
                               int mt )  -> std::pair< ParticleID, ParticleID > {

    unsigned int aa = std::round( AtomicMass( ma * neutronMass ).value );
    unsigned int ab = std::round( AtomicMass( mb * neutronMass ).value );

    std::string level = "_e0";
    if ( mt == 51 )  { level = "_e1"; }
    else if ( mt == 600 ) { level = "_e1"; };

    return { particleID.at( static_cast< unsigned int >( za ) * 1000 + aa ),
             elementID.at( static_cast< unsigned int >( zb ) )
               + std::to_string( ab )
               + level };
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
             residuals,
             endfPairs.Q() );
}
