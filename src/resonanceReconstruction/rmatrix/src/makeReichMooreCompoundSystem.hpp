CompoundSystem< ReichMoore, ShiftFactor >
makeReichMooreCompoundSystem(
    const ENDF::resolved::ReichMoore& endfReichMoore,
    const AtomicMass& neutronMass,
    const ElectricalCharge& elementaryCharge,
    const ParticleID& incident,
    const ParticleID& target,
    const std::optional< ChannelRadiusTable >& nro,
    unsigned int naps ) {

  // get some information
  unsigned int nls = endfReichMoore.numberLValues();
  unsigned int nlsc = endfReichMoore.numberLValuesForConvergence();
  double awri = endfReichMoore.lValues().front().atomicWeightRatio();
  double spin = endfReichMoore.spin();
  double ap = endfReichMoore.scatteringRadius();

  // particle pair
  ParticlePair in( Particle( incident, neutronMass,
                             0.0 * coulombs, 0.5, +1),
                   Particle( target, awri * neutronMass,
                             0.0 * coulombs, spin, +1) );

  // channel radii
  ChannelRadii radii = makeChannelRadii( ap, nro, naps,
                                         awri, neutronMass.value );

  // create the spin groups and assign the quantum numbers
  auto required = makeQuantumNumbers( in, nlsc );
  auto lvalues = endfReichMoore.lValues();
  std::vector< ParticleChannelData > data;
  for ( const auto& lvalue : lvalues ) {

    auto stuff = makeReichMooreChannelData( lvalue, neutronMass,
                                            elementaryCharge,
                                            incident, target, required,
                                            spin, ap, nro, naps );
    data.insert( data.end(), stuff.begin(), stuff.end() );
  }

  // add the remaining required elastic channels with no resonances
  for ( const auto& number : required ) {

    Channel< Neutron > elastic( in, in, 0.0 * electronVolt, number, radii );
    data.emplace_back( elastic,
                       std::vector< Energy >{},
                       std::vector< ReducedWidth >{} );
  }

  // return the resulting compound system
  return CompoundSystem< ReichMoore, ShiftFactor >(
             std::move( consolidateChannelData( data ) ) );
}
