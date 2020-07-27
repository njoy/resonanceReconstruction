/**
*  @brief Return the minimal energy grid derived from the resonance parameters
*         for the spin group
 *
 *  The minimal energy grid consists of all resonance energy values (only the
 *  ones that are positive) along with the resonance energies +/- half the
 *  total width.
 */
std::vector< Energy > grid() const {

  std::vector< Energy > grid;
  grid.reserve( 3 * this->resonanceTable().numberResonances() );

  // some useful lambdas
  auto toWidth = [] ( const auto& gamma, const auto& p ) -> Width
                    { return 2. * p * gamma * gamma; };
  auto totalWidth = [&] ( const auto& resonance ) -> Width {

    auto widths = ranges::view::zip_with(
                      toWidth,
                      resonance.widths(),
                      this->penetrabilities( resonance.energy() ) );

    Width total = 2. * resonance.eliminatedWidth()
                     * resonance.eliminatedWidth();
    return ranges::accumulate( widths, total );
  };

  for ( const auto& resonance : this->resonanceTable().resonances() ) {

    auto total = totalWidth( resonance );
    auto energy = resonance.energy();
    grid.push_back( energy - 0.5 * total );
    grid.push_back( energy );
    grid.push_back( energy + 0.5 * total );
  }
  grid |= ranges::action::sort | ranges::action::unique;

  return grid;
}
