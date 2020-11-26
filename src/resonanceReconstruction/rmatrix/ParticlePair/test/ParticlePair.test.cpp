#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction/rmatrix/ParticlePair.hpp"

using namespace njoy::resonanceReconstruction;
using namespace dimwits::constant;

// convenience typedefs
using Particle = rmatrix::Particle;
using ParticlePair = rmatrix::ParticlePair;
using Spin = rmatrix::Spin;
using Parity = rmatrix::Parity;
using ParticleID = rmatrix::ParticleID;
using ParticlePairID = rmatrix::ParticlePairID;

constexpr ElectricalCharge elementary = dimwits::constant::elementaryCharge;
constexpr double e = 1.6021766208e-19;

SCENARIO( "ParticlePair" ) {

  GIVEN( "valid data for a ParticlePair" ) {

    // particles
    Particle photon( ParticleID( "g" ), 0.0 * daltons, 0.0 * coulombs, 1., +1);
    Particle neutron( ParticleID( "n" ), 1.00866491582 * daltons,
                      0.0 * coulombs, 0.5, +1);
    Particle proton( ParticleID( "p" ), 1.00727647 * daltons,
                     elementary, 0.5, +1);
    Particle cl36( ParticleID( "Cl36_e0" ), 35.968306822 * daltons,
                   17.0 * elementary, 0., +1);
    Particle cl35( ParticleID( "Cl35_e0" ), 34.968852694 * daltons,
                   17.0 * elementary, 1.5, +1);
    Particle s36( ParticleID( "S36_e0" ), 35.967080699 * daltons,
                  16.0 * elementary, 1.5, +1);

    // custom particle pair identifier
    ParticlePairID id( "fission" );

    THEN( "a ParticlePair can be constructed without a pair ID" ) {

      ParticlePair pair( neutron, cl35 );

      CHECK( 1.00866491582 == Approx( pair.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair.particle().spin() ) );
      CHECK( +1 == pair.particle().parity() );

      CHECK( 34.968852694 == Approx( pair.residual().mass().value ) );
      CHECK( 17.0 * e == Approx( pair.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair.residual().spin() ) );
      CHECK( +1 == pair.residual().parity() );

      CHECK( "n,Cl35" == pair.pairID().symbol() );

      CHECK( 34.968852694 / ( 34.968852694 + 1.00866491582 )
             == Approx( pair.massRatio() ) );
      CHECK( 1.00866491582 * 34.968852694 / ( 34.968852694 + 1.00866491582 )
             == Approx( pair.reducedMass().value ) );
    } // THEN

    THEN( "a ParticlePair can be constructed with a pair ID" ) {

      ParticlePair pair( neutron, cl35, id );

      CHECK( 1.00866491582 == Approx( pair.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair.particle().spin() ) );
      CHECK( +1 == pair.particle().parity() );

      CHECK( 34.968852694 == Approx( pair.residual().mass().value ) );
      CHECK( 17.0 * e == Approx( pair.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair.residual().spin() ) );
      CHECK( +1 == pair.residual().parity() );

      CHECK( "fission" == pair.pairID().symbol() );

      CHECK( 34.968852694 / ( 34.968852694 + 1.00866491582 )
             == Approx( pair.massRatio() ) );
      CHECK( 1.00866491582 * 34.968852694 / ( 34.968852694 + 1.00866491582 )
             == Approx( pair.reducedMass().value ) );
    } // THEN
  } // GIVEN
} // SCENARIO
