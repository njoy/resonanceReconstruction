/**
*  @brief Return the minimal energy grid derived from the resonance parameters
*         for the spin group
 *
 *  The minimal energy grid consists of all resonance energy values (only the
 *  ones that are positive) along with the resonance energies +/- half the
 *  total width (including the eliminated width if it is defined).
 */
std::vector< Energy > grid() const {

  std::vector< Energy > grid;
  grid.reserve( 3 * this->resonanceTable().numberResonances() );

  // some useful lambdas
  auto toWidth = [] ( const auto& gamma, const auto& p ) -> Width {

    return 2. * p * gamma * gamma;
  };

  auto totalWidth = [&] ( const auto& resonance ) -> Width {

    auto widths = ranges::views::zip_with(
                      toWidth,
                      resonance.widths(),
                      this->penetrabilities( resonance.energy() ) );

    Width total = 2. * resonance.eliminatedWidth()
                     * resonance.eliminatedWidth();
    return ranges::accumulate( widths, total );
  };

  for ( const auto& resonance : this->resonanceTable().resonances() ) {

    decltype(auto) energy = resonance.energy();
    if ( energy > 0. * electronVolt ) {

      decltype(auto) total = totalWidth( resonance );
      grid.push_back( energy - 0.5 * total );
      grid.push_back( energy );
      grid.push_back( energy + 0.5 * total );
    }
  }

  ranges::cpp20::sort( grid );
  grid.erase( ranges::cpp20::unique( grid ), grid.end() );

  return grid;
}
