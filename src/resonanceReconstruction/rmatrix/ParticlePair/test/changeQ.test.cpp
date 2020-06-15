#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;
using namespace dimwits::constant;

// convenience typedefs
using Particle = rmatrix::Particle;
using ParticlePair = rmatrix::ParticlePair;
using Spin = rmatrix::Spin;
using Parity = rmatrix::Parity;
using ReactionID = rmatrix::ReactionID;

SCENARIO( "changeQ" ) {

  GIVEN( "valid data for a ParticlePair" ) {

    Particle neutron( "n", 1.008664 * daltons, 0.0 * coulombs, 0.5, +1 );
    Particle u235( "U235_e0", 235.0439299 * daltons,
                              92. * elementaryCharge, 0.0, +1 );
    QValue qValue = 0.0 * electronVolt;

    THEN( "the Q value can be changed" ) {

      ParticlePair pair( neutron, u235, qValue );
      REQUIRE( 0.0 == Approx( pair.Q().value ) );

      pair.changeQ( 10.0 * electronVolt );
      REQUIRE( 10.0 == Approx( pair.Q().value ) );
    } // THEN
  } // GIVEN
} // SCENARIO
