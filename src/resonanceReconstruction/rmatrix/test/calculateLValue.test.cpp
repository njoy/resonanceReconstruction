#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;
using namespace njoy::resonanceReconstruction::rmatrix;

SCENARIO( "calculateLValue" ) {

  GIVEN( "valid values for S, B and P" ) {

    THEN( "the appropriate value of L is calculated" ) {

      // SAMMY boundary condition
      REQUIRE( 0. == Approx( calculateLValue< ShiftFactor >( 1., 2., 3. ).real() ) );
      REQUIRE( 3. == Approx( calculateLValue< ShiftFactor >( 1., 2., 3. ).imag() ) );

      // Constant boundary condition
      REQUIRE( -1. == Approx( calculateLValue< Constant >( 1., 2., 3. ).real() ) );
      REQUIRE( 3. == Approx( calculateLValue< Constant >( 1., 2., 3. ).imag() ) );
    }
  } // GIVEN
} // SCENARIO

