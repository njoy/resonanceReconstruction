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

SCENARIO( "incident state" ) {

  GIVEN( "valid data for a ParticlePair" ) {

    Particle neutron( "n", 1.008664 * daltons, 0.0 * coulombs, 0.5, +1 );
    Particle u235( "U235_e0", 235.0439299 * daltons,
                              92. * elementaryCharge, 0.0, +1 );
    QValue qValue = 0.0 * electronVolt;

    THEN( "the incident state can be set at construction and toggled after" ) {

      ParticlePair pair( neutron, u235, qValue, true );
      REQUIRE( true == pair.incident() );

      pair.toggleIncident();
      REQUIRE( false == pair.incident() );

      pair.toggleIncident();
      REQUIRE( true == pair.incident() );
    } // THEN
  } // GIVEN
} // SCENARIO
