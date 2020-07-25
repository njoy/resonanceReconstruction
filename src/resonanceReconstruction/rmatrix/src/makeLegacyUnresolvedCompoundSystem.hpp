legacy::unresolved::CompoundSystem
makeLegacyUnresolvedCompoundSystem(
    const ENDF::unresolved::EnergyDependent& endfEnergyDependent,
    const AtomicMass& neutronMass,
    const ElectricalCharge& elementaryCharge,
    const ParticleID& incident,
    const ParticleID& target,
    unsigned int naps ) {

  // get some information
  unsigned int nls = endfEnergyDependent.NLS();
  double awri = endfEnergyDependent.lValues().front().atomicWeightRatio();
  double spin = endfEnergyDependent.spin();
  double ap = endfEnergyDependent.scatteringRadius();

  // particle pair
  ParticlePair in( Particle( incident, neutronMass,
                             0.0 * coulombs, 0.5, +1),
                   Particle( target, awri * neutronMass,
                             0.0 * coulombs, spin, +1) );

  // channel radii
  ChannelRadii radii = makeChannelRadii( ap, naps, awri, neutronMass.value );

  // some magic magic
  auto lvalues = endfEnergyDependent.lValues();
  auto repeatPair = ranges::view::repeat_n( in, nls );
  auto repeatRadii = ranges::view::repeat_n( radii, nls );
  std::vector< std::vector< legacy::unresolved::SpinGroup > > groups =
      ranges::view::zip_with(
          makeLegacyUnresolvedSpinGroups,
          lvalues,
          repeatPair,
          repeatRadii );

  return legacy::unresolved::CompoundSystem(
             std::move( groups | ranges::view::join ) );
}
