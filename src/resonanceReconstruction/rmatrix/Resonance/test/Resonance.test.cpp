#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using Resonance = rmatrix::Resonance;

SCENARIO( "Resonance" ) {

  GIVEN( "valid data for a Resonance" ) {

    // single resonance data
    Energy energy = 6.823616e+4 * electronVolt;
    ReducedWidth eliminated = 3.933600e-1 * rootElectronVolt;
    std::vector< ReducedWidth > widths = { 2.179040e+2 * rootElectronVolt,
                                           1.000000e-5 * rootElectronVolt };

    THEN( "a Resonance can be constructed" ) {
      Resonance resonance( energy, std::move( widths ), eliminated );

      REQUIRE( 6.823616e+4 == Approx( resonance.energy().value ) );

      REQUIRE( 3.933600e-1 == Approx( resonance.eliminatedWidth().value ) );

      REQUIRE( 2 == resonance.widths().size() );
      REQUIRE( 2.179040e+2 == Approx( resonance.widths()[0].value ) );
      REQUIRE( 1.000000e-5 == Approx( resonance.widths()[1].value ) );
    }
  } // GIVEN
} // SCENARIO


