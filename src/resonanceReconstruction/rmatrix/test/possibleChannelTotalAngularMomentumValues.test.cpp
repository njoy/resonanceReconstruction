#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;
using namespace njoy::resonanceReconstruction::rmatrix;

// convenience typedefs
using OrbitalAngularMomentum = rmatrix::OrbitalAngularMomentum;
using Spin = rmatrix::Spin;

SCENARIO( "getPossibleChannelTotalAngularMomentumValues" ) {

  GIVEN( "valid values for the orbital angular momentum, the particle and "
         "target spin" ) {

    THEN( "the appropriate total angular momentum values are generated" ) {

      // l=0, i=0.0, I=0.0, 0.5, 1.0
      auto values =
        rmatrix::getPossibleChannelTotalAngularMomentumValues( 0, 0.0, 0.0 );
      CHECK( 1 == Approx( values.size() ) );
      CHECK( 0.0 == Approx( values[0] ) );

      values =
        rmatrix::getPossibleChannelTotalAngularMomentumValues( 0, 0.0, 0.5 );
      CHECK( 1 == Approx( values.size() ) );
      CHECK( 0.5 == Approx( values[0] ) );

      values =
        rmatrix::getPossibleChannelTotalAngularMomentumValues( 0, 0.0, 1.0 );
      CHECK( 1 == Approx( values.size() ) );
      CHECK( 1.0 == Approx( values[0] ) );

      // l=0, i=0.5, I=0.0, 0.5, 1.0
      values =
        rmatrix::getPossibleChannelTotalAngularMomentumValues( 0, 0.5, 0.0 );
      CHECK( 1 == Approx( values.size() ) );
      CHECK( 0.5 == Approx( values[0] ) );

      values =
        rmatrix::getPossibleChannelTotalAngularMomentumValues( 0, 0.5, 0.5 );
      CHECK( 2 == Approx( values.size() ) );
      CHECK( 0.0 == Approx( values[0] ) );
      CHECK( 1.0 == Approx( values[1] ) );

      values =
        rmatrix::getPossibleChannelTotalAngularMomentumValues( 0, 0.5, 1.0 );
      CHECK( 2 == Approx( values.size() ) );
      CHECK( 0.5 == Approx( values[0] ) );
      CHECK( 1.5 == Approx( values[1] ) );

      // l=0, i=1.0, I=0.0, 0.5, 1.0
      values =
        rmatrix::getPossibleChannelTotalAngularMomentumValues( 0, 1.0, 0.0 );
      CHECK( 1 == Approx( values.size() ) );
      CHECK( 1.0 == Approx( values[0] ) );

      values =
        rmatrix::getPossibleChannelTotalAngularMomentumValues( 0, 1.0, 0.5 );
      CHECK( 2 == Approx( values.size() ) );
      CHECK( 0.5 == Approx( values[0] ) );
      CHECK( 1.5 == Approx( values[1] ) );

      values =
        rmatrix::getPossibleChannelTotalAngularMomentumValues( 0, 1.0, 1.0 );
      CHECK( 3 == Approx( values.size() ) );
      CHECK( 0.0 == Approx( values[0] ) );
      CHECK( 1.0 == Approx( values[1] ) );
      CHECK( 2.0 == Approx( values[2] ) );

      // l=1, i=0.0, I=0.0, 0.5, 1.0
      values =
        rmatrix::getPossibleChannelTotalAngularMomentumValues( 1, 0.0, 0.0 );
      CHECK( 1 == Approx( values.size() ) );
      CHECK( 1.0 == Approx( values[0] ) );

      values =
        rmatrix::getPossibleChannelTotalAngularMomentumValues( 1, 0.0, 0.5 );
      CHECK( 2 == Approx( values.size() ) );
      CHECK( 0.5 == Approx( values[0] ) );
      CHECK( 1.5 == Approx( values[1] ) );

      values =
        rmatrix::getPossibleChannelTotalAngularMomentumValues( 1, 0.0, 1.0 );
      CHECK( 3 == Approx( values.size() ) );
      CHECK( 0.0 == Approx( values[0] ) );
      CHECK( 1.0 == Approx( values[1] ) );
      CHECK( 2.0 == Approx( values[2] ) );

      // l=1, i=0.5, I=0.0, 0.5, 1.0
      values =
        rmatrix::getPossibleChannelTotalAngularMomentumValues( 1, 0.5, 0.0 );
      CHECK( 2 == Approx( values.size() ) );
      CHECK( 0.5 == Approx( values[0] ) );
      CHECK( 1.5 == Approx( values[1] ) );

      values =
        rmatrix::getPossibleChannelTotalAngularMomentumValues( 1, 0.5, 0.5 );
      CHECK( 3 == Approx( values.size() ) );
      CHECK( 0.0 == Approx( values[0] ) );
      CHECK( 1.0 == Approx( values[1] ) );
      CHECK( 2.0 == Approx( values[2] ) );

      values =
        rmatrix::getPossibleChannelTotalAngularMomentumValues( 1, 0.5, 1.0 );
      CHECK( 3 == Approx( values.size() ) );
      CHECK( 0.5 == Approx( values[0] ) );
      CHECK( 1.5 == Approx( values[1] ) );
      CHECK( 2.5 == Approx( values[2] ) );

      // l=1, i=1.0, I=0.0, 0.5, 1.0
      values =
        rmatrix::getPossibleChannelTotalAngularMomentumValues( 1, 1.0, 0.0 );
      CHECK( 3 == Approx( values.size() ) );
      CHECK( 0.0 == Approx( values[0] ) );
      CHECK( 1.0 == Approx( values[1] ) );
      CHECK( 2.0 == Approx( values[2] ) );

      values =
        rmatrix::getPossibleChannelTotalAngularMomentumValues( 1, 1.0, 0.5 );
      CHECK( 3 == Approx( values.size() ) );
      CHECK( 0.5 == Approx( values[0] ) );
      CHECK( 1.5 == Approx( values[1] ) );
      CHECK( 2.5 == Approx( values[2] ) );

      values =
        rmatrix::getPossibleChannelTotalAngularMomentumValues( 1, 1.0, 1.0 );
      CHECK( 3 == Approx( values.size() ) );
      CHECK( 1.0 == Approx( values[0] ) );
      CHECK( 2.0 == Approx( values[1] ) );
      CHECK( 3.0 == Approx( values[2] ) );
    }
  } // GIVEN
} // SCENARIO
