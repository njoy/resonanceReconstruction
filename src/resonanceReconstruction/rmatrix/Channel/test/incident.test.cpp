#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;
using namespace dimwits::constant;

// convenience typedefs
using Particle = rmatrix::Particle;
using ParticlePair = rmatrix::ParticlePair;
using Photon = rmatrix::Photon;
template < typename Type > using Channel = rmatrix::Channel< Type >;
using ChannelQuantumNumbers = rmatrix::ChannelQuantumNumbers;
using ChannelRadii = rmatrix::ChannelRadii;
using Photon = rmatrix::Photon;

SCENARIO( "incident state" ) {

  GIVEN( "valid data for a Channel" ) {

    Particle photon( "g", 0.0 * daltons, 0.0 * coulombs, 1., +1);
    Particle fe55( "Fe55_e0", 55.0 * daltons, 26.0 * coulombs, 0.0, +1);
    ParticlePair pair( photon, fe55, 0.0 * electronVolt, true );
    ChannelQuantumNumbers numbers( 0, 0.0, 0.5, +1 );
    ChannelRadii radii( 5.437300e-1 * rootBarn, 5.437300e-1 * rootBarn );

    THEN( "the incident state can be toggled" ) {

      Channel< Photon > channel( pair, numbers, radii, true );
      REQUIRE( true == channel.incident() );

      channel.toggleIncident();
      REQUIRE( false == channel.incident() );

      channel.toggleIncident();
      REQUIRE( true == channel.incident() );
    } // THEN
  } // GIVEN
} // SCENARIO
