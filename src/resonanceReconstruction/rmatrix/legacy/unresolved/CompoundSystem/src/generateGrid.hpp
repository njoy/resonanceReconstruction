/**
 *  @brief Generate the energy grid for which cross sections in the unresolved
 *         resonance region will be generated
 *
 *  As a first energy grid, this function will take every energy point for which
 *  unresolved resonance parameters are given. Normally, an evaluator is
 *  supposed to use the same energy grid for each l,J spin group. However, just
 *  in case data is used that does not conform this requirement, this function
 *  will merge all energy points found each l,J spin group and then generate a
 *  unique set of points.
 *
 *  @todo add additional points as required when the grid is too fine.
 */
static std::vector< Energy >
generateGrid( const std::vector< SpinGroup >& groups ) {

  // generate the initial grid using the energies of all resonance tables
  std::vector< Energy > grid  =
      groups | ranges::view::transform(
                   [] ( const auto& group )
                      { return group.resonanceTable().energies(); } )
             | ranges::view::join
             | ranges::to_vector
             | ranges::action::sort
             | ranges::action::unique;

  for ( const auto& energy : grid ) { std::cout << energy << std::endl; }

  return grid;
}
