legacy::unresolved::ResonanceTable
makeResonanceTable(
    const ENDF::unresolved::EnergyDependent::JValue& endfParameters ) {

  // usefull numbers

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
                           const ReducedWidth& capture,
                           const ReducedWidth& fission,
                           const ReducedWidth& competition ) {

    return Resonance( energy, spacing, elastic, capture, fission, competition );
  };

  // do some ranges magic
  std::vector< ChannelID > id = channels | ranges::view::transform( getID );
  if ( eliminated >= 0 ) {

    id.erase( id.begin() + eliminated );
  }
  auto resonances =
      ranges::view::zip_with( toResonance,
                              endfParameters.resonanceEnergies()
                                | ranges::view::transform( toEnergy ),
                              endfParameters.resonanceParameters() );

  return ResonanceTable( std::move( id ), std::move( resonances ) );
}
