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
    Width elastic = 1. * electronVolt;
    Width capture = 2. * electronVolt;
    Width fission = .5 * electronVolt;
    Width competition = 3. * electronVolt;
    double p = 2.;
    double q = 1.5;
    double s = -1.;

    THEN( "an Resonance can be constructed" ) {

      Resonance resonance( energy, elastic, capture, fission, competition,
                           p, q, s );

      CHECK( 1000. == Approx( resonance.energy().value ) );

      CHECK( 1. == Approx( resonance.elastic().value ) );
      CHECK( 2. == Approx( resonance.capture().value ) );
      CHECK( .5 == Approx( resonance.fission().value ) );
      CHECK( 3. == Approx( resonance.competition().value ) );

      CHECK( 1000.25 == Approx( resonance.energyPrime( -2. ).value ) );

      CHECK( .25 == Approx( resonance.elastic( 0.5 ).value ) );
      CHECK( 1. == Approx( resonance.competition( 0.5 ).value ) );
    } // THEN
  } // GIVEN
} // SCENARIO
