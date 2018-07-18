#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using Particle = rmatrix::Particle;
using ParticlePair = rmatrix::ParticlePair;
using Neutron = rmatrix::Neutron;
using Photon = rmatrix::Photon;
using Channel = rmatrix::Channel;
using ChannelID = rmatrix::ChannelID;
using ChannelQuantumNumbers = rmatrix::ChannelQuantumNumbers;
using ChannelRadii = rmatrix::ChannelRadii;
using ChannelType = rmatrix::ChannelType;
using BoundaryCondition = rmatrix::BoundaryCondition;
using Photon = rmatrix::Photon;

SCENARIO( "Channel" ) {

  GIVEN( "valid data for a Channel" ) {

    ChannelID id = "1";
    Particle photon( 0.0 * daltons, 0.0 * coulombs, 1., +1);
    Particle fe55( 55.0 * daltons, 26.0 * coulombs, 0.0, +1);
    ParticlePair pair( photon, fe55, 0.0 * electronVolt, "capture" );
    ChannelQuantumNumbers numbers( 0, 0.0, 0.5, +1 );
    ChannelRadii radii( 5.437300e-1 * rootBarn, 5.437300e-1 * rootBarn );
    BoundaryCondition boundary = 0.0;
    ChannelType type = Photon();

    THEN( "a Channel can be constructed" ) {

      Channel channel( id, pair, numbers, radii, boundary, type );

      REQUIRE( "1" == channel.channelID() );
      REQUIRE( 0.0 == Approx( channel.particlePair().particle().mass().value ) );
      REQUIRE( 0.0 == Approx( channel.particlePair().particle().charge().value ) );
      REQUIRE( 1.0 == Approx( channel.particlePair().particle().spin() ) );
      REQUIRE( +1 == channel.particlePair().particle().parity() );
      REQUIRE( 55. == Approx( channel.particlePair().residual().mass().value ) );
      REQUIRE( 26. == Approx( channel.particlePair().residual().charge().value ) );
      REQUIRE( 0.0 == Approx( channel.particlePair().residual().spin() ) );
      REQUIRE( +1 == channel.particlePair().residual().parity() );
      REQUIRE( 0.0 == Approx( channel.particlePair().Q().value ) );
      REQUIRE( "capture" == channel.particlePair().reaction() );
      REQUIRE( 0 == channel.quantumNumbers().orbitalAngularMomentum() );
      REQUIRE( 0.0 == channel.quantumNumbers().spin() );
      REQUIRE( 0.5 == channel.quantumNumbers().totalAngularMomentum() );
      REQUIRE( +1 == channel.quantumNumbers().parity() );
      REQUIRE( 5.437300e-1 ==
          Approx( channel.radii().penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      REQUIRE( 5.437300e-1 ==
          Approx( channel.radii().shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      REQUIRE( 5.437300e-1 ==
          Approx( channel.radii().phaseShiftRadius( 1e-5 * electronVolt ).value ) );
      REQUIRE( 0.0 == channel.boundaryCondition() );
    }
  } // GIVEN
} // SCENARIO


