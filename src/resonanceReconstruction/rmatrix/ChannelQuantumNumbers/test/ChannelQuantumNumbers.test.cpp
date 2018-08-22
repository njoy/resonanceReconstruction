#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using ChannelQuantumNumbers = rmatrix::ChannelQuantumNumbers;
using OrbitalAngularMomentum = rmatrix::OrbitalAngularMomentum;
using Spin = rmatrix::Spin;
using TotalAngularMomentum = rmatrix::TotalAngularMomentum;
using Parity = rmatrix::Parity;

SCENARIO( "ChannelQuantumNumbers" ) {

  GIVEN( "valid data for a ChannelQuantumNumbers" ) {

    OrbitalAngularMomentum l = 1;
    Spin s = 0.5;
    TotalAngularMomentum J = 1.5;
    Parity pi = +1;

    THEN( "a ChannelQuantumNumbers can be constructed" ) {
      ChannelQuantumNumbers numbers( l, s, J, pi );

      REQUIRE( 1 == numbers.orbitalAngularMomentum() );
      REQUIRE( 0.5 == numbers.spin() );
      REQUIRE( 1.5 == numbers.totalAngularMomentum() );
      REQUIRE( +1 == numbers.parity() );
      REQUIRE( "{1,1/2,3/2+}" == numbers.toString() );
    }
  } // GIVEN
} // SCENARIO


