SCENARIO("phaseShift"){
  auto awr = 2.360045E+2;

  auto a = channelRadius( awr );
  auto k = neutronWaveNumber( awr );
  auto phi = [&]( auto energy ){ return a(energy) * k(energy); };

  auto channelRatios =
    ranges::view::linear_distribute( -10., 30., 40 )
    | ranges::view::transform( []( double d ){ return std::pow( 2., d ); } )
    | ranges::view::transform( []( double d ){ return d * electronVolts; } )
    | ranges::view::transform( [&]( auto energy ){ return phi(energy); } )
    | ranges::to_vector;

  WHEN("l = 0"){
    auto l = 0;
    auto base = makeLvalue( l );

    const auto& lValue =
      static_cast< const breitWigner::singleLevel::Lvalue& >( base );

    auto trial =
      channelRatios
      | ranges::view::transform( [&]( auto ratio )
                                 { return lValue.phaseShift( ratio ); } );

    auto reference =
      channelRatios
      | ranges::view::transform( [&]( auto ratio )
                                 { return phaseShift( l, ratio ); } );

    for ( const auto& pair : ranges::view::zip( trial, reference ) ){
      auto trial = std::get<0>(pair);
      auto reference = std::get<1>(pair);
      REQUIRE( trial == reference );
    }
  }


  WHEN("l = 1"){
    auto l = 1;
    auto base = makeLvalue( l );

    const auto& lValue =
      static_cast< const breitWigner::singleLevel::Lvalue& >( base );

    auto trial =
      channelRatios
      | ranges::view::transform( [&]( auto ratio )
                                 { return lValue.phaseShift( ratio ); } );

    auto reference =
      channelRatios
      | ranges::view::transform( [&]( auto ratio )
                                 { return phaseShift( l, ratio ); } );

    for ( const auto& pair : ranges::view::zip( trial, reference ) ){
      auto trial = std::get<0>(pair);
      auto reference = std::get<1>(pair);
      REQUIRE( trial == reference );
    }
  }
}
