#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using Resonance = rmatrix::legacy::resolved::Resonance;

SCENARIO( "Resonance" ) {

  GIVEN( "valid data for a Resonance" ) {

    // single resonance data
    Energy energy = 1000. * electronVolt;
    ReducedWidth elastic = 1. * rootElectronVolt;
    ReducedWidth capture = 2. * rootElectronVolt;
    ReducedWidth fission = .5 * rootElectronVolt;
    ReducedWidth total = 5. * rootElectronVolt;

    THEN( "an Resonance can be constructed" ) {

      Resonance resonance( energy, total, elastic, capture, fission );

      CHECK( 1000. == Approx( resonance.energy().value ) );

      CHECK( 1. == Approx( resonance.elastic().value ) );
      CHECK( 2. == Approx( resonance.capture().value ) );
      CHECK( .5 == Approx( resonance.fission().value ) );
      CHECK( 5. == Approx( resonance.total().value ) );

      CHECK( .5 == Approx( resonance.elastic( 0.5 ).value ) );
      CHECK( 2. == Approx( resonance.competition( 0.5 ).value ) );
    } // THEN
  } // GIVEN
} // SCENARIO
