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

constexpr ElectricalCharge elementary = dimwits::constant::elementaryCharge;

SCENARIO( "ParticlePair" ) {

  GIVEN( "valid data for a ParticlePair" ) {

    // particles
    Particle photon( "g", 0.0 * daltons, 0.0 * coulombs, 1., +1);
    Particle neutron( "n", 1.00866491582 * daltons, 0.0 * coulombs, 0.5, +1);
    Particle proton( "p", 1.00727647 * daltons, elementary, 0.5, +1);
    Particle cl36( "Cl36_e0", 35.968306822 * daltons,
                              17.0 * elementary, 0., +1);
    Particle cl35( "Cl35_e0", 34.968852694 * daltons,
                              17.0 * elementary, 1.5, +1);
    Particle s36( "S36_e0", 35.967080699 * daltons,
                            16.0 * elementary, 1.5, +1);

    // custom particle pair identifier
    ParticlePairID id( "neutron" );

    THEN( "a ParticlePair can be constructed without a pair ID" ) {

      ParticlePair pair( neutron, cl35 );

      CHECK( 1.00866491582 == Approx( pair.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair.particle().spin() ) );
      CHECK( +1 == pair.particle().parity() );

      CHECK( 34.968852694 == Approx( pair.residual().mass().value ) );
      CHECK( 17.0 * 1.6021766208e-19
             == Approx( pair.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair.residual().spin() ) );
      CHECK( +1 == pair.residual().parity() );

      CHECK( "n,Cl35_e0" == pair.pairID() );

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
      CHECK( 17.0 * 1.6021766208e-19
             == Approx( pair.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair.residual().spin() ) );
      CHECK( +1 == pair.residual().parity() );

      CHECK( "neutron" == pair.pairID() );

      CHECK( 34.968852694 / ( 34.968852694 + 1.00866491582 )
             == Approx( pair.massRatio() ) );
      CHECK( 1.00866491582 * 34.968852694 / ( 34.968852694 + 1.00866491582 )
             == Approx( pair.reducedMass().value ) );
    } // THEN
  } // GIVEN
} // SCENARIO
