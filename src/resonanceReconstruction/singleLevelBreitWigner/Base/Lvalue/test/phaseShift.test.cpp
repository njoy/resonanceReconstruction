#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction::singleLevelBreitWigner;
using namespace dimwits;

Base::Lvalue makeLvalue();

SCENARIO("phaseShift"){
  auto lValue = makeLvalue();
  auto awr = 2.360045E+2;
  auto l = 0;
  
  auto a = channelRadius( awr );
  auto k = neutronWaveNumber( awr );
  auto phi = [&]( auto energy ){ return a(energy) * k(energy); };
  
  auto channelRatios =
    ranges::view::linear_distribute( -10., 30., 40 )
    | ranges::view::transform( []( double d ){ return std::pow( 2., d ); } )
    | ranges::view::transform( []( double d ){ return d * electronVolts; } )
    | ranges::view::transform( [&]( auto energy ){ return phi(energy); } )
    | ranges::to_vector;
  {
    auto trial =
      channelRatios
      | ranges::view::transform( [&]( auto ratio )
                                 { return lValue.phaseShift( ratio ); } );

    auto reference =
      channelRatios
      | ranges::view::transform( [&]( auto ratio )
                                 { return phaseShift( l, ratio ); } );

    RANGES_FOR( auto pair, ranges::view::zip( trial, reference ) ){
      auto trial = std::get<0>(pair);
      auto reference = std::get<1>(pair);
      REQUIRE( trial == reference );
    }
  }

  lValue.orbitalAngularMomentum = 1;
  l = 1;

  {
    auto trial =
      channelRatios
      | ranges::view::transform( [&]( auto ratio )
                                 { return lValue.phaseShift( ratio ); } );

    auto reference =
      channelRatios
      | ranges::view::transform( [&]( auto ratio )
                                 { return phaseShift( l, ratio ); } );

    RANGES_FOR( auto pair, ranges::view::zip( trial, reference ) ){
      auto trial = std::get<0>(pair);
      auto reference = std::get<1>(pair);
      REQUIRE( trial == reference );
    }
  }
}

