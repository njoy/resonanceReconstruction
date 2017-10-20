#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;
using namespace dimwits;

njoy::ENDFtk::section::Type<2> resonances();

namespace {

struct Apply : breitWigner::multilevel::Apply {
  using breitWigner::multilevel::Apply::build;
};

}

template< typename T >
void ignore( T&& ){}

SCENARIO( "Integration test" ){
  const auto Na22 = resonances();
  const auto& isotope = Na22.isotopes.front();
  const auto& energyRange = isotope.energyRanges().front();
  const auto& mlbw = std::experimental::get< 2 >( energyRange );
  const auto type = Apply().build( mlbw, channelRadius( 22. ), radius( mlbw.AP() ) );
  {
    const auto xs = type( 11129.000000001111 * electronVolts );
    REQUIRE( xs.capture.value  == Approx( 7.6715238015690504E-005 ) );
    REQUIRE( xs.elastic.value  == Approx( 4.5011595167555498 ) );
  }
  {
    const auto xs = type( 3.0000000000000001E-005 * electronVolts );
    REQUIRE( xs.capture.value  == Approx( 8479.8299771152415 ) );
    REQUIRE( xs.elastic.value  == Approx( 83.300829493115927 ) );
  }
}

std::string Sodium22Resonances();

njoy::ENDFtk::section::Type<2>
resonances(){
  auto string = Sodium22Resonances();
  auto start = string.begin();
  auto it = start;
  auto end = string.end();
  long lineNumber = 1;

  auto head = njoy::ENDFtk::HEAD( it, end, lineNumber );
  return njoy::ENDFtk::section::Type<2>( head, it, end, lineNumber, 1122 );
}

std::string Sodium22Resonances(){
  return njoy::utility::slurpFileToMemory( "Na-22-mt151" );
}
