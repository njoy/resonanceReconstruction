static std::vector< Energy >
makeEnergies( const std::vector< Resonance >& resonances ) {

  return ranges::to< std::vector< Energy > >(
         resonances
         | ranges::cpp20::views::transform(
               [] ( const auto& resonance ) -> decltype(auto)
                  { return resonance.energy(); } ) );
}

static std::vector< LevelSpacing >
makeSpacings( const std::vector< Resonance >& resonances ) {

  return ranges::to< std::vector< LevelSpacing > >(
         resonances
         | ranges::cpp20::views::transform(
               [] ( const auto& resonance ) -> decltype(auto)
                  { return resonance.levelSpacing(); } ) );
}

static std::vector< ReducedWidth >
makeElasticWidths( const std::vector< Resonance >& resonances ) {

  return ranges::to< std::vector< ReducedWidth > >(
         resonances
         | ranges::cpp20::views::transform(
               [] ( const auto& resonance ) -> decltype(auto)
                  { return resonance.elastic(); } ) );
}

static std::vector< Width >
makeCaptureWidths( const std::vector< Resonance >& resonances ) {

  return ranges::to< std::vector< Width > >(
         resonances
         | ranges::cpp20::views::transform(
               [] ( const auto& resonance ) -> decltype(auto)
                  { return resonance.capture(); } ) );
}

static std::vector< Width >
makeFissionWidths( const std::vector< Resonance >& resonances ) {

  return ranges::to< std::vector< Width > >(
         resonances
         | ranges::cpp20::views::transform(
               [] ( const auto& resonance ) -> decltype(auto)
                  { return resonance.fission(); } ) );
}

static std::vector< Width >
makeCompetitionWidths( const std::vector< Resonance >& resonances ) {

  return ranges::to< std::vector< Width > >(
         resonances
         | ranges::cpp20::views::transform(
               [] ( const auto& resonance ) -> decltype(auto)
                  { return resonance.competition(); } ) );
}

static std::tuple< LevelSpacingTable, ReducedWidthTable,
                   WidthTable, WidthTable, WidthTable >
makeTables( const std::vector< Resonance >& resonances ) {

  std::vector< Energy > energies = makeEnergies( resonances );
  std::vector< LevelSpacing > spacings = makeSpacings( resonances );
  std::vector< ReducedWidth > elastic = makeElasticWidths( resonances );
  std::vector< Width > capture = makeCaptureWidths( resonances );
  std::vector< Width > fission = makeFissionWidths( resonances );
  std::vector< Width > competition = makeCompetitionWidths( resonances );

  return { LevelSpacingTable( std::vector< Energy >( energies ),
                              std::move( spacings ) ),
           ReducedWidthTable( std::vector< Energy >( energies ),
                              std::move( elastic ) ),
           WidthTable( std::vector< Energy >( energies ),
                       std::move( capture ) ),
           WidthTable( std::vector< Energy >( energies ),
                       std::move( fission ) ),
           WidthTable( std::move( energies ), std::move( competition ) ) };
}
