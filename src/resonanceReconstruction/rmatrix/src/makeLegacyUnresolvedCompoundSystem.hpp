template < typename Unresolved >
std::vector< double > getEnergies( const Unresolved&,
                                   double lower, double upper ) {

  return { lower, upper };
}

template <> std::vector< double >
getEnergies(
  const endf::UnresolvedEnergyDependentFissionWidths& unresolved,
  double, double ) {

  return ranges::to< std::vector< double > >( unresolved.energies() );
}

template < typename Unresolved >
int getUnresolvedInterpolation( const Unresolved& ) { return 5; }

template <>
int getUnresolvedInterpolation(
        const endf::UnresolvedEnergyDependent& unresolved ) {

  auto verifyInterpolation = [] ( const auto& interpolants ) {

    int interpolation = interpolants[0];
    if ( ranges::cpp20::count( interpolants, interpolation )
           != ranges::cpp20::distance( interpolants ) ) {

      throw std::runtime_error( "Different interpolation schemes for "
                                "unresolved l,J values are currently not "
                                "implemented" );
    }
    return interpolation;
  };

  auto getInterpolation = [&] ( const auto& values ) {

    return verifyInterpolation(
               values | ranges::views::transform(
                            [] ( const auto& jvalue )
                               { return jvalue.INT(); } ) );
  };

  return verifyInterpolation(
             unresolved.lValues()
                 | ranges::views::transform(
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
  auto repeatPair = ranges::views::repeat_n( in, nls );
  auto repeatRadii = ranges::views::repeat_n( radii, nls );
  auto repeatEnergies = ranges::views::repeat_n(
                            getEnergies( unresolved, lower, upper ), nls );

  std::vector< legacy::unresolved::SpinGroup > groups;
  auto data = ranges::views::zip_with(
                  makeLegacyUnresolvedSpinGroups< typename Unresolved::LValue >,
                  lvalues,
                  repeatPair,
                  repeatRadii,
                  repeatEnergies );
  for ( const auto& entry : data ) {

    groups.insert( groups.end(), entry.begin(), entry.end() );
  }

  return legacy::unresolved::CompoundSystem( std::move( groups ), interpolation );
}
