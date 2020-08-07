template < typename Unresolved >
std::vector< double > getEnergies( const Unresolved& unresolved,
                                   double lower, double upper ) {

  return { lower, upper };
}

template <> std::vector< double >
getEnergies(
  const ENDF::unresolved::EnergyDependentFissionWidths& unresolved,
  double lower, double upper ) {

  return unresolved.energies();
}

template < typename Unresolved >
legacy::unresolved::CompoundSystem
makeLegacyUnresolvedCompoundSystem(
    const Unresolved& unresolved,
    const AtomicMass& neutronMass,
    const ElectricalCharge& elementaryCharge,
    const ParticleID& incident,
    const ParticleID& target,
    unsigned int naps,
    double lower,
    double upper ) {

  // get some information
  unsigned int nls = unresolved.NLS();
  double awri = unresolved.lValues().front().atomicWeightRatio();
  double spin = unresolved.spin();
  double ap = unresolved.scatteringRadius();

  // particle pair
  ParticlePair in( Particle( incident, neutronMass,
                             0.0 * coulombs, 0.5, +1),
                   Particle( target, awri * neutronMass,
                             0.0 * coulombs, spin, +1) );

  // channel radii
  ChannelRadii radii = makeChannelRadii( ap, naps, awri, neutronMass.value );

  // some magic magic
  auto lvalues = unresolved.lValues();
  auto repeatPair = ranges::view::repeat_n( in, nls );
  auto repeatRadii = ranges::view::repeat_n( radii, nls );
  auto repeatEnergies = ranges::view::repeat_n(
                            getEnergies( unresolved, lower, upper ), nls );
  std::vector< std::vector< legacy::unresolved::SpinGroup > > groups =
      ranges::view::zip_with(
          makeLegacyUnresolvedSpinGroups< typename Unresolved::LValue >,
          lvalues,
          repeatPair,
          repeatRadii,
          repeatEnergies );

  return legacy::unresolved::CompoundSystem(
             std::move( groups | ranges::view::join ) );
}
