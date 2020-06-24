#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using ChannelRadii = rmatrix::ChannelRadii;

SCENARIO( "ChannelRadii" ) {

  GIVEN( "valid data for a ChannelRadii" ) {

    ChannelRadius radius = 0.1 * rootBarns;

    ChannelRadius trueRadius = 0.2 * rootBarns;
    ChannelRadius effectiveRadius = 0.3 * rootBarns;

    ChannelRadius pRadius = 0.4 * rootBarns;
    ChannelRadius sRadius = 0.5 * rootBarns;
    ChannelRadius phiRadius = 0.6 * rootBarns;

    THEN( "a ChannelRadii can be constructed using constant radii" ) {

      Energy energy = 1e-5 * electronVolts;

      ChannelRadii radii1( radius );

      CHECK( 0.1 == Approx( radii1.penetrabilityRadius( energy ).value ) );
      CHECK( 0.1 == Approx( radii1.shiftFactorRadius( energy ).value ) );
      CHECK( 0.1 == Approx( radii1.phaseShiftRadius( energy ).value ) );

      ChannelRadii radii2( trueRadius, effectiveRadius );

      CHECK( 0.2 == Approx( radii2.penetrabilityRadius( energy ).value ) );
      CHECK( 0.2 == Approx( radii2.shiftFactorRadius( energy ).value ) );
      CHECK( 0.3 == Approx( radii2.phaseShiftRadius( energy ).value ) );

      ChannelRadii radii3( pRadius, sRadius, phiRadius );

      CHECK( 0.4 == Approx( radii3.penetrabilityRadius( energy ).value ) );
      CHECK( 0.5 == Approx( radii3.shiftFactorRadius( energy ).value ) );
      CHECK( 0.6 == Approx( radii3.phaseShiftRadius( energy ).value ) );
    } // THEN
  } // GIVEN
} // SCENARIO
