#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;
using namespace dimwits;

njoy::ENDFtk::section::Type<2> resonances();

namespace {

struct Apply : reichMoore::Apply {
  using reichMoore::Apply::build;
};

}

template< typename T >
void ignore( T&& ){}

SCENARIO( "Integration test" ){
  const auto Al27 = resonances();
  const auto& isotope = Al27.isotopes.front();
  const auto& energyRange = isotope.energyRanges().front();
  const auto& rm = std::experimental::get< 3 >( energyRange );
  njoy::Log::info(
    "\n Al-27"
    "\n============="
    "\n LRU: {}"
    "\n LRF: {}"
    "\n NRO: {}"
    "\n NAPS: {}"
    "\n AP: {}\n\n",
    rm.LRU(), rm.LRF(), rm.NRO(), rm.NAPS(), rm.AP() );

  const auto type =
    Apply().build( rm, reichMoore::Both{}, radius( rm.AP() ), true );

  {
    const auto xs = type( 1.E-5 * electronVolts );
    REQUIRE( xs.elastic.value == Approx( 1.4250586145968933 ) );
    REQUIRE( xs.fission.value == Approx( 0.0000000000000000 ) );
    REQUIRE( xs.capture.value == Approx( 11.742759183335906 ) );
  }{
    const auto xs = type( 1.E-4 * electronVolts );
    REQUIRE( xs.elastic.value == Approx( 1.4250586080911192 ) );
    REQUIRE( xs.fission.value == Approx( 0.0000000000000000 ) );
    REQUIRE( xs.capture.value == Approx( 3.7133864261889338 ) );
  }{
    const auto xs = type( 1.E-3 * electronVolts );
    REQUIRE( xs.elastic.value == Approx( 1.4250585448655608 ) );
    REQUIRE( xs.fission.value == Approx( 0.0000000000000000 ) );
    REQUIRE( xs.capture.value == Approx( 1.1742756516065123 ) );
  }{
    const auto xs = type( 844927.20000008447 * electronVolts );
    REQUIRE( xs.elastic.value == Approx( 4.0296898644697023 ) );
    REQUIRE( xs.fission.value == Approx( 0.0000000000000000 ) );
    REQUIRE( xs.capture.value == Approx( 1.6742377842945329E-004 ) );
  }
}

njoy::ENDFtk::section::Type<2>
resonances(){
  auto endfFile = njoy::utility::slurpFileToMemory( "n-013_Al_027.endf" );

  auto begin = endfFile.begin();
  auto end = endfFile.end();

  njoy::ENDFtk::syntaxTree::Tape< std::string::iterator > tape( begin, end );

  auto& material = *( tape.begin() );

  auto MAT = material.MAT();
  long lineNumber = 1;
  return material
         .fileNumber(2)
         .sectionNumber(151).parse<2>( lineNumber, MAT );
}
