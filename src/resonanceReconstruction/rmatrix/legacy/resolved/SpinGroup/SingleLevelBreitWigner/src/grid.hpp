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

  for ( const auto& resonance : this->resonanceTable().resonances() ) {

    auto energy = resonance.energy();
    if ( energy > 0. * electronVolt ) {

      auto total = resonance.total();
      grid.push_back( energy - 0.5 * total );
      grid.push_back( energy );
      grid.push_back( energy + 0.5 * total );
    }
  }
  grid |= ranges::action::sort | ranges::action::unique;

  return grid;
}
