template < typename Formalism, typename BoundaryOption >
inline SpinGroup< Formalism, BoundaryOption >
makeSpinGroup(
    const ParticlePair& incident,
    const std::vector< ParticlePair > pairs,
    const ENDF::resolved::RMatrixLimited::ParticlePairs& endfPairs,
    const ENDF::resolved::RMatrixLimited::SpinGroup& endfSpinGroup,
    bool reducedWidthsFlag,
    Formalism,
    BoundaryOption ) {

  std::vector< ParticleChannel > channels =
      makeParticleChannels( incident, pairs, endfPairs,
                            endfSpinGroup.channels() );
  ResonanceTable table = makeResonanceTable( channels,
                                             endfSpinGroup.parameters(),
                                             reducedWidthsFlag );

  return SpinGroup< Formalism, BoundaryOption >( std::move( channels ),
                                                 std::move( table ) );
}

template < typename BoundaryOption >
inline SpinGroup< ReichMoore, BoundaryOption >
makeSpinGroup(
    const ParticlePair& incident,
    const std::vector< ParticlePair > pairs,
    const ENDF::resolved::RMatrixLimited::ParticlePairs& endfPairs,
    const ENDF::resolved::RMatrixLimited::SpinGroup& endfSpinGroup,
    bool reducedWidthsFlag,
    ReichMoore,
    BoundaryOption ) {

  // get channels and determine the eliminated channel index
  std::vector< ParticleChannel > channels =
      makeParticleChannels( incident, pairs, endfPairs,
                            endfSpinGroup.channels() );
  unsigned int eliminated = [&] {

    auto getParticlePair = [] ( const ParticleChannel& channel ) {

      return std::visit( [] ( const auto& channel )
                            { return channel.particlePair(); },
                         channel );
    };

    auto eliminatedPair = pairs[ rmatrix::eliminated( endfPairs ) ];
    auto pairs = channels | ranges::view::transform( getParticlePair );
    auto found = std::find_if(
                     ranges::begin( pairs ),
                     ranges::end( pairs ),
                     [&] ( const auto& pair )
                         { return eliminatedPair.pairID() == pair.pairID(); } );
    return std::distance( ranges::begin( pairs ), found );
  }();

  ResonanceTable table = makeResonanceTable( channels,
                                             endfSpinGroup.parameters(),
                                             reducedWidthsFlag,
                                             eliminated );
  channels.erase( channels.begin() + eliminated );
  return SpinGroup< ReichMoore, BoundaryOption >( std::move( channels ),
                                                  std::move( table ) );
}
