template < typename Formalism, typename BoundaryOption >
inline CompoundSystem< Formalism, BoundaryOption >
makeCompoundSystem(
    const ENDF::resolved::RMatrixLimited& endfRMatrix,
    const AtomicMass& neutronMass,
    const ElectricalCharge& elementaryCharge,
    Formalism formalism,
    BoundaryOption boundaryOption ) {

  auto endfPairs = endfRMatrix.particlePairs();
  auto pairs = makeParticlePairs( endfPairs, neutronMass, elementaryCharge );
  auto incident = pairs[ rmatrix::incident( endfPairs ) ];
  bool reducedWidthsFlag = endfRMatrix.reducedWidths() == 0 ? false : true;

  std::vector< SpinGroup< Formalism, BoundaryOption > > spingroups =
      endfRMatrix.spinGroups()
        | ranges::view::transform(
              [&] ( const auto& spingroup ) {

                return makeSpinGroup( incident, pairs,
                                      endfRMatrix.particlePairs(),
                                      spingroup, reducedWidthsFlag,
                                      formalism, boundaryOption );
              } );

  return CompoundSystem< Formalism, BoundaryOption >( std::move( spingroups ) );
}
