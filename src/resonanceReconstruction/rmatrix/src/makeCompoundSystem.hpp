template < typename Formalism, typename BoundaryOption >
CompoundSystem< Formalism, BoundaryOption >
makeCompoundSystem(
    const endf::RMatrixLimited& endfRMatrix,
    const AtomicMass& neutronMass,
    const ElectricalCharge& elementaryCharge,
    const ParticleID& incident,
    const ParticleID& target,
    Formalism formalism,
    BoundaryOption boundaryOption ) {

  auto endfPairs = endfRMatrix.particlePairs();
  auto pairs = makeParticlePairs( endfPairs, neutronMass, elementaryCharge,
                                  incident, target );
  auto in = pairs[ rmatrix::incident( endfPairs ) ];
  bool reducedWidthsFlag = endfRMatrix.reducedWidths() == 0 ? false : true;
  auto spingroups = endfRMatrix.spinGroups();

  auto data = spingroups
    | ranges::view::transform(
          [&] ( const auto& spingroup )
              { return makeParticleChannelData(
                           in, pairs, endfPairs, spingroup, reducedWidthsFlag,
                           formalism, boundaryOption ); } )
    | ranges::to_vector;
  std::vector< ParticleChannelData > channels =
      consolidateChannelData( data | ranges::view::join );

  return CompoundSystem< Formalism, BoundaryOption >( std::move( channels ) );
}
