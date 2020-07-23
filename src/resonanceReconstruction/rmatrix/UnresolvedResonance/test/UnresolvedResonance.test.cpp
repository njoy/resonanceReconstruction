#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using UnresolvedResonance = rmatrix::UnresolvedResonance;

SCENARIO( "UnresolvedResonance" ) {

  GIVEN( "valid data for a UnresolvedResonance" ) {

    // single resonance data
    Energy energy = 6.823616e+4 * electronVolt;
    Energy spacing = 3.933600e-1 * electronVolt;
    std::vector< ReducedWidth > widths = { 2.179040e+2 * rootElectronVolt,
                                           1.000000e-5 * rootElectronVolt };

    THEN( "an UnresolvedResonance can be constructed" ) {

      UnresolvedResonance resonance( energy, spacing, std::move( widths ) );

      CHECK( 6.823616e+4 == Approx( resonance.energy().value ) );

      CHECK( 3.933600e-1 == Approx( resonance.levelSpacing().value ) );

      CHECK( 2 == resonance.widths().size() );
      CHECK( 2.179040e+2 == Approx( resonance.widths()[0].value ) );
      CHECK( 1.000000e-5 == Approx( resonance.widths()[1].value ) );
    } // THEN
  } // GIVEN
} // SCENARIO
