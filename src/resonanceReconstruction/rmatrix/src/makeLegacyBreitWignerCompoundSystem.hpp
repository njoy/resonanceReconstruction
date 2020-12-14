template < typename Formalism, typename BreitWigner >
legacy::resolved::CompoundSystem< Formalism >
makeLegacyBreitWignerCompoundSystem(
    const BreitWigner& endfBreitWigner,
    const AtomicMass& neutronMass,
    const ParticleID& incident,
    const ParticleID& target,
    const std::optional< ChannelRadiusTable >& nro,
    unsigned int naps,
    Formalism formalism ) {

  // get some information
  unsigned int nls = endfBreitWigner.numberLValues();
  double awri = endfBreitWigner.lValues().front().atomicWeightRatio();
  double spin = endfBreitWigner.spin();
  double ap = endfBreitWigner.scatteringRadius();

  // particle pair
  ParticlePair in( Particle( incident, neutronMass,
                             0.0 * coulombs, 0.5, +1),
                   Particle( target,
                             awri * neutronMass,
                             0.0 * coulombs, spin, +1) );

  // create the spin groups and assign the quantum numbers
  auto required = makeQuantumNumbers( in, nls );
  auto lvalues = endfBreitWigner.lValues();
  std::vector< legacy::resolved::SpinGroup< Formalism > > groups;
  for ( const auto& lvalue : lvalues ) {

    auto stuff = makeLegacyBreitWignerSpinGroups(
                     lvalue, neutronMass, incident, target,
                     required, spin, ap, nro, naps, formalism );
    groups.insert( groups.end(), stuff.begin(), stuff.end() );
  }

  // return the resulting compound system
  return legacy::resolved::CompoundSystem< Formalism >( std::move( groups ) );
}
