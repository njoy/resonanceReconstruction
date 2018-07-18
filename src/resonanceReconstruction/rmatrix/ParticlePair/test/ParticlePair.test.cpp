#define CATCH_CONFIG_MAIN

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

SCENARIO( "ParticlePair" ) {

  GIVEN( "valid data for a ParticlePair" ) {

    // particle
    Particle neutron( 1.008664 * daltons, 0.0 * coulombs, 0.5, +1 );
    Particle u235( 235.0439299 * daltons, 92. * elementaryCharge, 0.0, +1 );
    QValue qValue = 0.0 * electronVolt;
    ReactionID reaction( "elastic" );

    THEN( "a ParticlePair can be constructed" ) {
      ParticlePair pair( neutron, u235, qValue, reaction );

      REQUIRE( 1.008664 == Approx( pair.particle().mass().value ) );
      REQUIRE( 0.0 == Approx( pair.particle().charge().value ) );
      REQUIRE( 0.5 == Approx( pair.particle().spin() ) );
      REQUIRE( +1 == pair.particle().parity() );

      REQUIRE( 235.0439299 == Approx( pair.residual().mass().value ) );
      REQUIRE( 0.0 == Approx( pair.residual().charge().value ) );
      REQUIRE( 0.0 == Approx( pair.residual().spin() ) );
      REQUIRE( +1 == pair.residual().parity() );

      REQUIRE( 0.0 == Approx( pair.Q().value ) );
      REQUIRE( "elastic" == pair.reaction() );

      REQUIRE( 1.00435393058671 == Approx( pair.reducedMass().value ) );
      REQUIRE( 6.912508e-06 == Approx(
        pair.waveNumber( 1e-5 * electronVolt ).value ) );
      REQUIRE( 0.0 == Approx( pair.etaParameter( 1e-5 * electronVolt ).value ) );
    }
  } // GIVEN
} // SCENARIO

