/**
 *  @brief Generate the energy grid for which cross sections in the unresolved
 *         resonance region will be generated
 *
 *  As a first energy grid, this function will take every energy point for which
 *  unresolved resonance parameters are given. Normally, an evaluator is
 *  supposed to use the same energy grid for each l,J spin group. However, just
 *  in case data is used that does not conform this requirement, this function
 *  will merge all energy points found in each l,J spin group and then generate
 *  a unique set of points.
 *
 *  After that, the function will go over consectutive points and verify if
 *  the points vary by more than 25%. If they do, points will be inserted using
 *  a fixed 13 point per decade scheme.
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

  // go over each energy point (no need to do this for empty grids or grids that
  // have only 1 element in them because an exception will be thrown anyway)
  if ( grid.size() > 1 ) {

    static constexpr std::array< double, 13 > points =
      { 1.0, 1.25, 1.5, 1.7, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.2, 8.5 };

    auto iter = grid.begin();
    auto previous = iter++;
    for ( ; iter != grid.end(); ++iter ) {

      if ( *iter / *previous > 1.25 ) {

        int exponent = std::floor( std::log10( std::abs( ( *previous ).value ) ) );
        if ( points.back() * std::pow( 10., exponent ) * electronVolt <= *previous ) {

          ++exponent;
        }

        auto values =
            points | ranges::view::transform(
                         [&] ( const auto& point ) -> Energy
                             { return point * std::pow( 10., exponent )
                                      * electronVolt; } );

        auto begin = std::upper_bound( ranges::begin( values ), ranges::end( values ), *previous );
        auto end = std::lower_bound( ranges::begin( values ), ranges::end( values ), *iter );
        auto distance = std::distance( begin, end );
        iter = grid.insert( iter, begin, end );
        std::advance( iter, distance - 1 );
      }
      previous = iter;
    }
  }

  for ( const auto& energy : grid ) { std::cout << energy << std::endl; }

  return grid;
}
