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
using ParticlePairID = rmatrix::ParticlePairID;

SCENARIO( "ParticlePair" ) {

  GIVEN( "valid data for a ParticlePair" ) {

    // particle
    Particle neutron( "n", 1.008664 * daltons, 0.0 * coulombs, 0.5, +1 );
    Particle u235( "U235_e0", 235.0439299 * daltons,
                              92. * elementaryCharge, 0.0, +1 );
    QValue qValue = 0.0 * electronVolt;
    ParticlePairID id( "fission" );

    THEN( "a ParticlePair can be constructed without a pair ID" ) {

      ParticlePair pair( neutron, u235, qValue );

      REQUIRE( 1.008664 == Approx( pair.particle().mass().value ) );
      REQUIRE( 0.0 == Approx( pair.particle().charge().value ) );
      REQUIRE( 0.5 == Approx( pair.particle().spin() ) );
      REQUIRE( +1 == pair.particle().parity() );

      REQUIRE( 235.0439299 == Approx( pair.residual().mass().value ) );
      REQUIRE( 0.0 == Approx( pair.residual().charge().value ) );
      REQUIRE( 0.0 == Approx( pair.residual().spin() ) );
      REQUIRE( +1 == pair.residual().parity() );

      REQUIRE( 235.0439299 / ( 235.0439299 + 1.008664 )
               == Approx( pair.massRatio() ) );

      REQUIRE( 0.0 == Approx( pair.Q().value ) );
      REQUIRE( "n,U235_e0" == pair.pairID() );

      REQUIRE( false == pair.incident() );

      REQUIRE( 1.00435393058671 == Approx( pair.reducedMass().value ) );
      REQUIRE( 6.9172282659633E-06 == Approx(
        pair.waveNumber( 1e-5 * electronVolt ).value ) );
      REQUIRE( 0.0 == Approx(
        pair.sommerfeldParameter( 1e-5 * electronVolt ) ) );
    } // THEN

    THEN( "a ParticlePair can be constructed with a pair ID" ) {

      ParticlePair pair( neutron, u235, qValue, id );

      REQUIRE( 1.008664 == Approx( pair.particle().mass().value ) );
      REQUIRE( 0.0 == Approx( pair.particle().charge().value ) );
      REQUIRE( 0.5 == Approx( pair.particle().spin() ) );
      REQUIRE( +1 == pair.particle().parity() );

      REQUIRE( 235.0439299 == Approx( pair.residual().mass().value ) );
      REQUIRE( 0.0 == Approx( pair.residual().charge().value ) );
      REQUIRE( 0.0 == Approx( pair.residual().spin() ) );
      REQUIRE( +1 == pair.residual().parity() );

      REQUIRE( 235.0439299 / ( 235.0439299 + 1.008664 )
               == Approx( pair.massRatio() ) );

      REQUIRE( 0.0 == Approx( pair.Q().value ) );
      REQUIRE( "fission" == pair.pairID() );

      REQUIRE( false == pair.incident() );

      REQUIRE( 1.00435393058671 == Approx( pair.reducedMass().value ) );
      REQUIRE( 6.9172282659633E-06 == Approx(
        pair.waveNumber( 1e-5 * electronVolt ).value ) );
      REQUIRE( 0.0 == Approx(
        pair.sommerfeldParameter( 1e-5 * electronVolt ) ) );
    } // THEN
  } // GIVEN
} // SCENARIO
