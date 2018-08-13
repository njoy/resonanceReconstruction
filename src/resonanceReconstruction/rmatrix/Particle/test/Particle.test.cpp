#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using Particle = rmatrix::Particle;
using Spin = rmatrix::Spin;
using Parity = rmatrix::Parity;

SCENARIO( "Particle" ) {

  GIVEN( "valid data for a Particle" ) {

    // neutron
    AtomicMass neutronMass = 1.008664 * daltons;
    ElectricalCharge neutronCharge = 0.0 * coulombs;
    Spin neutronSpin = 0.5;
    Parity neutronParity = +1;

    // proton
    AtomicMass protonMass = 1.007276 * daltons;
    ElectricalCharge protonCharge = dimwits::constant::elementaryCharge;
    Spin protonSpin = 0.5;
    Parity protonParity = +1;

    // U235
    AtomicMass u235Mass = 235.0439299 * daltons;
    ElectricalCharge u235Charge = 0.0 * coulombs;
    Spin u235Spin = 0.0;
    Parity u235Parity = +1;

    THEN( "a Particle can be constructed" ) {
      Particle neutron( "n", neutronMass, neutronCharge,
                             neutronSpin, neutronParity );

      REQUIRE( "n" == neutron.particleID() );
      REQUIRE( 1.008664 == Approx( neutron.mass().value ) );
      REQUIRE( 0.0 == Approx( neutron.charge().value ) );
      REQUIRE( 0.5 == Approx( neutron.spin() ) );
      REQUIRE( +1 == neutron.parity() );

      Particle proton( "p", protonMass, protonCharge,
                            protonSpin, protonParity );

      REQUIRE( "p" == proton.particleID() );
      REQUIRE( 1.007276 == Approx( proton.mass().value ) );
      REQUIRE( 1.60217662e-19 == Approx( proton.charge().value ) );
      REQUIRE( 0.5 == Approx( proton.spin() ) );
      REQUIRE( +1 == proton.parity() );

      Particle u235( "U235", u235Mass, u235Charge, u235Spin, u235Parity );

      REQUIRE( "U235" == u235.particleID() );
      REQUIRE( 235.0439299 == Approx( u235.mass().value ) );
      REQUIRE( 0.0 == Approx( u235.charge().value ) );
      REQUIRE( 0.0 == Approx( u235.spin() ) );
      REQUIRE( +1 == u235.parity() );
    }
  } // GIVEN
} // SCENARIO


