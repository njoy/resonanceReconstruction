legacy::unresolved::ResonanceTable
makeLegacyUnresolvedResonanceTable(
    const endf::UnresolvedEnergyDependent::JValue& endfParameters,
    const std::vector< double >& ) {

  // some usefull lambdas
  auto toEnergy = [&] ( double value ) -> Energy {

    return value * electronVolt;
  };
  auto toWidth = [&] ( double value ) -> Width {

    return value * electronVolt;
  };
  auto toReducedWidth = [&] ( double value ) -> ReducedWidth {

    return value * rootElectronVolt;
  };
  auto toLevelSpacing = [&] ( double value ) -> LevelSpacing {

    return value * electronVolt;
  };
  auto toResonance = [&] ( const Energy& energy,
                           const LevelSpacing& spacing,
                           const ReducedWidth& elastic,
                           const Width& capture,
                           const Width& fission,
                           const Width& competition ) {

    return legacy::unresolved::Resonance( energy, spacing, elastic,
                                          capture, fission, competition );
  };

  // do some ranges magic
  std::vector< legacy::unresolved::Resonance > resonances =
    ranges::to< std::vector< legacy::unresolved::Resonance > >(
      ranges::views::zip_with( toResonance,
                              endfParameters.energies()
                                | ranges::views::transform( toEnergy ),
                              endfParameters.averageLevelSpacings()
                                | ranges::views::transform( toLevelSpacing ),
                              endfParameters.averageNeutronWidths()
                                | ranges::views::transform( toReducedWidth ),
                              endfParameters.averageGammaWidths()
                                | ranges::views::transform( toWidth ),
                              endfParameters.averageFissionWidths()
                                | ranges::views::transform( toWidth ),
                              endfParameters.averageCompetitiveWidths()
                                | ranges::views::transform( toWidth ) ) );

  return legacy::unresolved::ResonanceTable(
             std::move( resonances ),
             { static_cast<unsigned int>( endfParameters.neutronWidthDegreesFreedom() ),
               static_cast<unsigned int>( endfParameters.gammaWidthDegreesFreedom() ),
               static_cast<unsigned int>( endfParameters.fissionWidthDegreesFreedom() ),
               static_cast<unsigned int>( endfParameters.competitiveWidthDegreesFreedom() ) } );
}

template< typename Range >
legacy::unresolved::ResonanceTable
makeLegacyUnresolvedResonanceTable(
    const endf::UnresolvedEnergyIndependent::JValue< Range >& endfParameters,
    const std::vector< double >& energies ) {

  // some usefull lambdas
  auto toEnergy = [&] ( double value ) -> Energy {

    return value * electronVolt;
  };
  auto toWidth = [&] ( double value ) -> Width {

    return value * electronVolt;
  };
  auto toReducedWidth = [&] ( double value ) -> ReducedWidth {

    return value * rootElectronVolt;
  };
  auto toLevelSpacing = [&] ( double value ) -> LevelSpacing {

    return value * electronVolt;
  };
  auto toResonance = [&] ( const Energy& energy,
                           const LevelSpacing& spacing,
                           const ReducedWidth& elastic,
                           const Width& capture,
                           const Width& fission,
                           const Width& competition ) {

    return legacy::unresolved::Resonance( energy, spacing, elastic,
                                          capture, fission, competition );
  };

  double lower = energies.front();
  double upper = energies.back();
  auto spacing = toLevelSpacing( endfParameters.averageLevelSpacing() );
  auto elastic = toReducedWidth( endfParameters.averageNeutronWidth() );
  auto gamma = toWidth( endfParameters.averageGammaWidth() );
  auto fission = toWidth( endfParameters.averageFissionWidth() );
  auto competition = toWidth( endfParameters.averageCompetitiveWidth() );
  return legacy::unresolved::ResonanceTable(
             { toResonance( toEnergy( lower ), spacing, elastic,
                            gamma, fission, competition ),
               toResonance( toEnergy( upper ), spacing, elastic,
                            gamma, fission, competition ) },
             { static_cast<unsigned int>( endfParameters.neutronWidthDegreesFreedom() ),
               static_cast<unsigned int>( endfParameters.gammaWidthDegreesFreedom() ),
               static_cast<unsigned int>( endfParameters.fissionWidthDegreesFreedom() ),
               static_cast<unsigned int>( endfParameters.competitiveWidthDegreesFreedom() ) } );
}

legacy::unresolved::ResonanceTable
makeLegacyUnresolvedResonanceTable(
    const endf::UnresolvedEnergyDependentFissionWidths::JValue& endfParameters,
    const std::vector< double >& energies ) {

  // some usefull lambdas
  auto toEnergy = [&] ( double value ) -> Energy {

    return value * electronVolt;
  };
  auto toWidth = [&] ( double value ) -> Width {

    return value * electronVolt;
  };
  auto toReducedWidth = [&] ( double value ) -> ReducedWidth {

    return value * rootElectronVolt;
  };
  auto toLevelSpacing = [&] ( double value ) -> LevelSpacing {

    return value * electronVolt;
  };
  auto toResonance = [&] ( const Energy& energy,
                           const LevelSpacing& spacing,
                           const ReducedWidth& elastic,
                           const Width& capture,
                           const Width& fission,
                           const Width& competition ) {

    return legacy::unresolved::Resonance( energy, spacing, elastic,
                                          capture, fission, competition );
  };

  // do some ranges magic
  auto ne = energies.size();
  auto spacing = toLevelSpacing( endfParameters.averageLevelSpacing() );
  auto elastic = toReducedWidth( endfParameters.averageNeutronWidth() );
  auto gamma = toWidth( endfParameters.averageGammaWidth() );
  auto competition = toWidth( endfParameters.averageCompetitiveWidth() );
  std::vector< legacy::unresolved::Resonance > resonances =
    ranges::to< std::vector< legacy::unresolved::Resonance > >(
      ranges::views::zip_with( toResonance,
                              energies | ranges::views::transform( toEnergy ),
                              ranges::views::repeat_n( spacing, ne ),
                              ranges::views::repeat_n( elastic, ne ),
                              ranges::views::repeat_n( gamma, ne ),
                              endfParameters.averageFissionWidths()
                                | ranges::views::transform( toWidth ),
                              ranges::views::repeat_n( competition, ne ) ) );

  return legacy::unresolved::ResonanceTable(
             std::move( resonances ),
             { static_cast<unsigned int>( endfParameters.neutronWidthDegreesFreedom() ),
               static_cast<unsigned int>( endfParameters.gammaWidthDegreesFreedom() ),
               static_cast<unsigned int>( endfParameters.fissionWidthDegreesFreedom() ),
               static_cast<unsigned int>( endfParameters.competitiveWidthDegreesFreedom() ) } );
}
