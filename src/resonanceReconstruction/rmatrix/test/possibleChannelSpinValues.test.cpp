#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;
using namespace njoy::resonanceReconstruction::rmatrix;

// convenience typedefs
using Spin = rmatrix::Spin;

SCENARIO( "possibleChannelSpinValues" ) {

  GIVEN( "valid values for the particle and target spin" ) {

    THEN( "the appropriate channel spin values are generated" ) {

      // i=0.0, I=0.0, 0.5, 1.0
      auto values = rmatrix::possibleChannelSpinValues( 0.0, 0.0 );
      REQUIRE( 1 == Approx( values.size() ) );
      REQUIRE( 0.0 == Approx( values[0] ) );

      values = rmatrix::possibleChannelSpinValues( 0.0, 0.5 );
      REQUIRE( 1 == Approx( values.size() ) );
      REQUIRE( 0.5 == Approx( values[0] ) );

      values = rmatrix::possibleChannelSpinValues( 0.0, 1.0 );
      REQUIRE( 1 == Approx( values.size() ) );
      REQUIRE( 1.0 == Approx( values[0] ) );

      // i=0.5, I=0.0, 0.5, 1.0
      values = rmatrix::possibleChannelSpinValues( 0.5, 0.0 );
      REQUIRE( 1 == Approx( values.size() ) );
      REQUIRE( 0.5 == Approx( values[0] ) );

      values = rmatrix::possibleChannelSpinValues( 0.5, 0.5 );
      REQUIRE( 2 == Approx( values.size() ) );
      REQUIRE( 0.0 == Approx( values[0] ) );
      REQUIRE( 1.0 == Approx( values[1] ) );

      values = rmatrix::possibleChannelSpinValues( 0.5, 1.0 );
      REQUIRE( 2 == Approx( values.size() ) );
      REQUIRE( 0.5 == Approx( values[0] ) );
      REQUIRE( 1.5 == Approx( values[1] ) );

      // i=1.0, I=0.0, 0.5, 1.0
      values = rmatrix::possibleChannelSpinValues( 1.0, 0.0 );
      REQUIRE( 1 == Approx( values.size() ) );
      REQUIRE( 1.0 == Approx( values[0] ) );

      values = rmatrix::possibleChannelSpinValues( 1.0, 0.5 );
      REQUIRE( 2 == Approx( values.size() ) );
      REQUIRE( 0.5 == Approx( values[0] ) );
      REQUIRE( 1.5 == Approx( values[1] ) );

      values = rmatrix::possibleChannelSpinValues( 1.0, 1.0 );
      REQUIRE( 3 == Approx( values.size() ) );
      REQUIRE( 0.0 == Approx( values[0] ) );
      REQUIRE( 1.0 == Approx( values[1] ) );
      REQUIRE( 2.0 == Approx( values[2] ) );
    }
  } // GIVEN
} // SCENARIO

