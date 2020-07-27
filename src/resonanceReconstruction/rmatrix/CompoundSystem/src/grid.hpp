/**
 *  @brief Return the minimal energy grid derived from the resonance parameters
 *         for the entire compound system
 *
 *  The minimal energy grid consists of all resonance energy values (only the
 *  ones that are positive) along with the resonance energies +/- half the
 *  total width.
 */
std::vector< Energy > grid() const {

  // generate the minimal grid using the minimal grid from each spin group
  std::vector< Energy > grid;
  for ( const auto& group : this->spinGroups() ) {

    auto groupgrid = group.grid();
    grid.insert( grid.end(), groupgrid.begin(), groupgrid.end() );
  }
  grid |= ranges::action::sort | ranges::action::unique;

  return grid;
}
