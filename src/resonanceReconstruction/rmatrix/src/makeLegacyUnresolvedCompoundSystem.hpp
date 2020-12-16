template < typename Unresolved >
std::vector< double > getEnergies( const Unresolved&,
                                   double lower, double upper ) {

  return { lower, upper };
}

template <> std::vector< double >
getEnergies(
  const endf::UnresolvedEnergyDependentFissionWidths& unresolved,
  double, double ) {

  return unresolved.energies();
}

template < typename Unresolved >
int getUnresolvedInterpolation( const Unresolved& ) { return 5; }

template <>
int getUnresolvedInterpolation(
        const endf::UnresolvedEnergyDependent& unresolved ) {

  auto verifyInterpolation = [] ( const auto& interpolants ) {

    int interpolation = interpolants[0];
    if ( ranges::count( interpolants, interpolation )
           != ranges::distance( interpolants ) ) {

      throw std::runtime_error( "Different interpolation schemes for "
                                "unresolved l,J values are currently not "
                                "implemented" );
    }
    return interpolation;
  };

  auto getInterpolation = [&] ( const auto& values ) {

    return verifyInterpolation(
               values | ranges::view::transform(
                            [] ( const auto& jvalue )
                               { return jvalue.INT(); } ) );
  };

  return verifyInterpolation(
             unresolved.lValues()
                 | ranges::view::transform(
                       [&] ( const auto& lvalue )
                           { return getInterpolation( lvalue.jValues() ); } ) );
}

template < typename Unresolved >
legacy::unresolved::CompoundSystem
makeLegacyUnresolvedCompoundSystem(
    const Unresolved& unresolved,
    const AtomicMass& neutronMass,
    const ParticleID& incident,
    const ParticleID& target,
    const std::optional< ChannelRadiusTable >& nro,
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
  ChannelRadii radii = makeChannelRadii( ap, nro, naps,
                                         awri, neutronMass.value );

  // interpolation parameter
  auto interpolation = getUnresolvedInterpolation( unresolved );

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
             std::move( groups | ranges::view::join ), interpolation );
}
