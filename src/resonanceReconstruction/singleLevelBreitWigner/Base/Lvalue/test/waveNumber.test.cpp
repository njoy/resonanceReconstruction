#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction::singleLevelBreitWigner;
using namespace dimwits;

Base::Lvalue makeLvalue();

SCENARIO("wave number"){
  auto lValue = makeLvalue();

  auto energies =
    ranges::view::linear_distribute( -10., 30., 40 )
    | ranges::view::transform( []( double d ){ return std::pow( 2., d ); } )
    | ranges::view::transform( []( double d ){ return d * electronVolts; } );

  RANGES_FOR( auto energy, energies ){
    auto trial = lValue.waveNumber( energy );
    auto reference = neutronWaveNumber( 2.360045E+2 )( energy );
    REQUIRE( trial == reference );
  }
}

