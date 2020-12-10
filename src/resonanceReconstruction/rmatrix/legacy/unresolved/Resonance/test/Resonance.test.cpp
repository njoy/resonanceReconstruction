#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction/rmatrix/legacy/unresolved/Resonance.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using Resonance = rmatrix::legacy::unresolved::Resonance;

SCENARIO( "Resonance" ) {

  GIVEN( "valid data for a Resonance" ) {

    // single resonance data
    Energy energy = 6.823616e+4 * electronVolt;
    LevelSpacing spacing = 3.933600e-1 * electronVolt;
    ReducedWidth elastic = 2.179040e+2 * rootElectronVolt;
    Width capture = 1.000000e-5 * electronVolt;
    Width fission = 0. * electronVolt;
    Width competition = 0.5 * electronVolt;

    THEN( "an Resonance can be constructed" ) {

      Resonance resonance( energy, spacing, elastic,
                           capture, fission, competition );

      CHECK( 6.823616e+4 == Approx( resonance.energy().value ) );

      CHECK( 3.933600e-1 == Approx( resonance.levelSpacing().value ) );

      CHECK( 2.179040e+2 == Approx( resonance.elastic().value ) );
      CHECK( 1.000000e-5 == Approx( resonance.capture().value ) );
      CHECK( 0. == Approx( resonance.fission().value ) );
      CHECK( 0.5 == Approx( resonance.competition().value ) );
    } // THEN
  } // GIVEN
} // SCENARIO
