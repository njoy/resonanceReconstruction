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

      CHECK( "n" == neutron.particleID() );
      CHECK( 1.008664 == Approx( neutron.mass().value ) );
      CHECK( 0.0 == Approx( neutron.charge().value ) );
      CHECK( 0.5 == Approx( neutron.spin() ) );
      CHECK( +1 == neutron.parity() );

      Particle proton( "p", protonMass, protonCharge,
                            protonSpin, protonParity );

      CHECK( "p" == proton.particleID() );
      CHECK( 1.007276 == Approx( proton.mass().value ) );
      CHECK( 1.60217662e-19 == Approx( proton.charge().value ) );
      CHECK( 0.5 == Approx( proton.spin() ) );
      CHECK( +1 == proton.parity() );

      Particle u235( "U235", u235Mass, u235Charge, u235Spin, u235Parity );

      CHECK( "U235" == u235.particleID() );
      CHECK( 235.0439299 == Approx( u235.mass().value ) );
      CHECK( 0.0 == Approx( u235.charge().value ) );
      CHECK( 0.0 == Approx( u235.spin() ) );
      CHECK( +1 == u235.parity() );
    } // THEN
  } // GIVEN

  GIVEN( "invalid data for a Particle" ) {

    // neutron
    AtomicMass neutronMass = 1.008664 * daltons;
    AtomicMass negativeMass = -1.008664 * daltons;
    ElectricalCharge neutronCharge = 0.0 * coulombs;
    ElectricalCharge negativeCharge = -1.0 * coulombs;
    Spin neutronSpin = 0.5;
    Parity neutronParity = +1;

    THEN( "an exception is thrown" ) {

      // negative mass
      CHECK_THROWS( Particle( "n", negativeMass, neutronCharge,
                              neutronSpin, neutronParity ) );

      // negative charge
      CHECK_THROWS( Particle( "n", neutronMass, negativeCharge,
                              neutronSpin, neutronParity ) );
    } // THEN
  } // GIVEN
} // SCENARIO
