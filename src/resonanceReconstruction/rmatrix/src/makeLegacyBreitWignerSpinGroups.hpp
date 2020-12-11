template < typename Formalism >
std::vector< legacy::resolved::SpinGroup< Formalism > >
makeLegacyBreitWignerSpinGroups(
    const endf::BreitWignerLValue& endfLValue,
    const AtomicMass& neutronMass,
    const ParticleID& incident,
    const ParticleID& target,
    std::vector< ChannelQuantumNumbers >& available,
    double spin,
    double ap,
    const std::optional< ChannelRadiusTable >& nro,
    unsigned int naps,
    Formalism ) {

  // the spin groups
  std::vector< legacy::resolved::SpinGroup< Formalism > > groups;

  // useful numbers
  unsigned int l = endfLValue.orbitalMomentum();
  auto awri = endfLValue.atomicWeightRatio();

  // useful lambdas
  auto toResonance = [] ( const auto& energy, const auto& elastic,
                          const auto& capture, const auto& fission,
                          const auto& competition, const auto& p,
                          const auto& q, const auto& s ) {

    return legacy::resolved::Resonance( energy, elastic, capture, fission,
                                        competition, p, q, s );
  };

  // particle pair
  ParticlePair in( Particle( incident, neutronMass,
                             0.0 * coulombs, 0.5, +1),
                   Particle( target,
                             awri * neutronMass,
                             0.0 * coulombs, spin, +1) );

  // channel radii
  ChannelRadii radii = makeChannelRadii( ap, nro, naps,
                                         awri, neutronMass.value );

  // the different J values referenced in the table
  std::vector< double > jvalues = endfLValue.spinValues()
                                      | ranges::to_vector
                                      | ranges::action::sort
                                      | ranges::action::unique;

  // loop over each J value
  auto resonances = endfLValue.resonances();
  for ( double J : jvalues ) {

    // filter out resonances with this J value
    auto current = resonances | ranges::view::filter(
                                     [J] ( const auto& resonance )
                                         { return resonance.spin() == J; } );

    // create the elastic channel
    Channel< Neutron > cElastic( in, in, 0. * electronVolt,
                                 retrieveQuantumNumber( l, J, available ),
                                 radii );

    // calculate P, Q and S
    unsigned int nr = ranges::distance( current );
    auto qx = endfLValue.QX() * electronVolt;
    bool lrx = endfLValue.competitiveWidthFlag();
    std::vector< Energy > energies =
        current | ranges::view::transform(
                    [] ( const auto& resonance )
                       { return resonance.resonanceEnergy() * electronVolt; } );
    std::vector< double >  p =
        energies | ranges::view::transform(
                     [&] ( const auto& energy )
                         { return cElastic.penetrability( energy ); } );
    std::vector< double > q =
        lrx ? energies
                | ranges::view::transform(
                    [&] ( const auto& energy )
                        { return cElastic.penetrability( energy - qx ); } )
            : p;
    std::vector< double > s =
        energies | ranges::view::transform(
                     [&] ( const auto& energy )
                         { return cElastic.shiftFactor( energy ); } );

    // get the widths
    auto elastic =
          current | ranges::view::transform(
                      [] ( const auto& resonance )
                         { return resonance.neutronWidth() * electronVolt; } );
    auto capture =
          current | ranges::view::transform(
                      [] ( const auto& resonance )
                         { return resonance.gammaWidth() * electronVolt; } );
    auto fission =
          current | ranges::view::transform(
                      [] ( const auto& resonance )
                         { return resonance.fissionWidth() * electronVolt; } );
    std::vector< Width > competition =
    lrx ? current
              | ranges::view::transform(
                  [&] ( const auto& resonance )
                      { return resonance.competitiveWidth() * electronVolt; } )
        : std::vector< Width >( nr, 0. * electronVolt );

    // make the spin group
    groups.emplace_back( std::move( cElastic ),
                         legacy::resolved::ResonanceTable(
                             ranges::view::zip_with(
                                 toResonance,
                                 energies, elastic, capture, fission,
                                 competition, p, q, s ) ),
                         qx );
  }

  // return the groups
  return groups;
}
