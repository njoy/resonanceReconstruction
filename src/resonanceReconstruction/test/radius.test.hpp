SCENARIO("radius"){
  SECTION("double overload"){
    const double AP = 10.0;

    auto reference = ranges::views::repeat_n( AP * rootBarns, 10 );

    auto trial = ranges::views::linear_distribute( -10, 30, 40 )
      | ranges::views::transform( []( double d ){ return std::pow( 2.0, d ); } )
      | ranges::views::transform( []( double d ){ return d * electronVolts ; } )
      | ranges::views::transform
        ( [ ap = radius( AP ) ]( auto&& e ){ return ap( e ); } );

    for ( const auto pair : ranges::views::zip( trial, reference ) ){
      const auto trial = std::get<0>(pair);
      const auto reference = std::get<1>(pair);
      REQUIRE( trial == reference );
    }
  }

  SECTION("TAB1 overload"){
    auto tab1 = []{
      std::vector< long > boundaries = { 3, 6 };
      std::vector< long > interpolants = { 1, 2 };
      std::vector< double > energies = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
      std::vector< double > radii = { 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 };

      return endf::ScatteringRadius( std::move( boundaries ),
                                     std::move( interpolants ),
                                     std::move( energies ),
                                     std::move( radii ) );
    }();

    auto trial = ranges::views::linear_distribute( 1.0, 6.0, 11 )
      | ranges::views::transform( []( double e ){ return e * electronVolts; } )
      | ranges::views::transform
        ( [ ape = radius( tab1 ) ]( auto&& e ){ return ape( e ); } );

    std::vector< Quantity<RootBarn> > reference = { 3.0 * rootBarn,
                                                    3.0 * rootBarn,
                                                    4.0 * rootBarn,
                                                    4.0 * rootBarn,
                                                    5.0 * rootBarn,
                                                    5.5 * rootBarn,
                                                    6.0 * rootBarn,
                                                    6.5 * rootBarn,
                                                    7.0 * rootBarn,
                                                    7.5 * rootBarn,
                                                    8.0 * rootBarn };

    for ( const auto pair : ranges::views::zip( trial, reference ) ){
      const auto trial = std::get<0>(pair);
      const auto reference = std::get<1>(pair);
      REQUIRE( trial.value == Approx(reference.value) );
    }
  }
}
