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

SCENARIO( "special functions" ) {

  GIVEN( "Two ParticlePair with/without charged particles" ) {

    // particle
    Particle neutron( "n", 1.008664 * daltons, 0.0 * coulombs, 0.5, +1 );
    Particle proton( "p", 1.007276 * daltons, elementaryCharge, 0.5, +1 );
    Particle u235( "U235_e0", 235.0439299 * daltons,
                              92. * elementaryCharge, 0.0, +1 );
    QValue qValue = 0.0 * electronVolt;

    ParticlePair neutronPair( neutron, u235, qValue );
    ParticlePair protonPair( proton, u235, qValue );

    THEN( "the special function give the right answers" ) {

      REQUIRE( 1.00435393058671 == Approx( neutronPair.reducedMass().value ) );
      REQUIRE( 6.917228e-06 == Approx(
        neutronPair.waveNumber( 1e-5 * electronVolt ).value ) );
      REQUIRE( 0.0 == Approx(
        neutronPair.etaParameter( 1e-5 * electronVolt ).value ) );

      REQUIRE( 1.00297775913185 == Approx( protonPair.reducedMass().value ) );
      REQUIRE( 6.917228e-06 == Approx(
        protonPair.waveNumber( 1e-5 * electronVolt ).value ) );
      REQUIRE( 5.395568e-38 == Approx(
        protonPair.etaParameter( 1e-5 * electronVolt ).value ) );
    }
  } // GIVEN
} // SCENARIO

