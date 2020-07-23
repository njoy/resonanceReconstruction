static std::vector< Energy >
makeEnergies( const std::vector< UnresolvedResonance >& resonances ) {

  return resonances
         | ranges::view::transform( [] ( const auto& resonance )
                                       { return resonance.energy(); } );
}

static std::vector< Energy >
makeSpacings( const std::vector< UnresolvedResonance >& resonances ) {

  return resonances
         | ranges::view::transform( [] ( const auto& resonance )
                                       { return resonance.levelSpacing(); } );
}

static std::vector< ReducedWidth >
makeWidths( unsigned int i,
            const std::vector< UnresolvedResonance >& resonances ) {

  return resonances
         | ranges::view::transform( [i] ( const auto& resonance )
                                        { return resonance.widths()[i]; } );
}

static Table< Energy, Energy >
makeLevelSpacingTable( const std::vector< UnresolvedResonance >& resonances ) {

  std::vector< Energy > energies = makeEnergies( resonances );
  std::vector< Energy > spacings = makeSpacings( resonances );
  return Table< Energy, Energy >( std::move( energies ), std::move( spacings ) );
}

static std::vector< Table< Energy, ReducedWidth > >
makeWidthTables( const std::vector< UnresolvedResonance >& resonances ) {

  std::vector< Table< Energy, ReducedWidth > > tables;
  unsigned int size = resonances.front().widths().size();
  for ( unsigned int i = 0; i < size; ++i ) {

    std::vector< Energy > energies = makeEnergies( resonances );
    std::vector< ReducedWidth > widths = makeWidths( i, resonances );
    tables.emplace_back( std::move( energies ), std::move( widths ) );
  }
  return tables;
}
