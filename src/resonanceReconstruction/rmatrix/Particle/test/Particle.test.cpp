#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction/rmatrix/Particle.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using Particle = rmatrix::Particle;
using ParticleID = rmatrix::ParticleID;
using Spin = rmatrix::Spin;
using Parity = rmatrix::Parity;

SCENARIO( "Particle" ) {

  GIVEN( "valid data for a Particle" ) {

    // neutron
    ParticleID neutronID( "n" );
    AtomicMass neutronMass = 1.008664 * daltons;
    ElectricalCharge neutronCharge = 0.0 * coulombs;
    Spin neutronSpin = 0.5;
    Parity neutronParity = +1;

    // proton
    ParticleID protonID( "p" );
    AtomicMass protonMass = 1.007276 * daltons;
    ElectricalCharge protonCharge = dimwits::constant::elementaryCharge;
    Spin protonSpin = 0.5;
    Parity protonParity = +1;

    // U235
    ParticleID u235ID( "U235" );
    AtomicMass u235Mass = 235.0439299 * daltons;
    ElectricalCharge u235Charge = 92 * dimwits::constant::elementaryCharge;
    Spin u235Spin = 0.0;
    Parity u235Parity = +1;

    THEN( "a Particle can be constructed" ) {

      Particle neutron( neutronID, neutronMass, neutronCharge,
                        neutronSpin, neutronParity );

      CHECK( "n" == neutron.particleID().symbol() );
      CHECK( 1.008664 == Approx( neutron.mass().value ) );
      CHECK( 0.0 == Approx( neutron.charge().value ) );
      CHECK( 0.5 == Approx( neutron.spin() ) );
      CHECK( +1 == neutron.parity() );

      neutron = Particle( neutronID, neutronMass,
                          neutronSpin, neutronParity );

      CHECK( "n" == neutron.particleID().symbol() );
      CHECK( 1.008664 == Approx( neutron.mass().value ) );
      CHECK( 0.0 == Approx( neutron.charge().value ) );
      CHECK( 0.5 == Approx( neutron.spin() ) );
      CHECK( +1 == neutron.parity() );

      Particle proton( protonID, protonMass, protonCharge,
                       protonSpin, protonParity );

      CHECK( "p" == proton.particleID().symbol() );
      CHECK( 1.007276 == Approx( proton.mass().value ) );
      CHECK( 1.60217662e-19 == Approx( proton.charge().value ) );
      CHECK( 0.5 == Approx( proton.spin() ) );
      CHECK( +1 == proton.parity() );

      proton = Particle( protonID, protonMass,
                         protonSpin, protonParity );

      CHECK( "p" == proton.particleID().symbol() );
      CHECK( 1.007276 == Approx( proton.mass().value ) );
      CHECK( 1.60217662e-19 == Approx( proton.charge().value ) );
      CHECK( 0.5 == Approx( proton.spin() ) );
      CHECK( +1 == proton.parity() );

      Particle u235( u235ID, u235Mass, u235Charge, u235Spin, u235Parity );

      CHECK( "U235" == u235.particleID().symbol() );
      CHECK( 235.0439299 == Approx( u235.mass().value ) );
      CHECK( 92 * 1.60217662e-19 == Approx( u235.charge().value ) );
      CHECK( 0.0 == Approx( u235.spin() ) );
      CHECK( +1 == u235.parity() );

      u235 = Particle( u235ID, u235Mass, u235Spin, u235Parity );

      CHECK( "U235" == u235.particleID().symbol() );
      CHECK( 235.0439299 == Approx( u235.mass().value ) );
      CHECK( 92 * 1.60217662e-19 == Approx( u235.charge().value ) );
      CHECK( 0.0 == Approx( u235.spin() ) );
      CHECK( +1 == u235.parity() );
    } // THEN
  } // GIVEN

  GIVEN( "invalid data for a Particle" ) {

    // neutron
    ParticleID neutronID( "n" );
    AtomicMass neutronMass = 1.008664 * daltons;
    AtomicMass negativeMass = -1.008664 * daltons;
    ElectricalCharge neutronCharge = 0.0 * coulombs;
    ElectricalCharge negativeCharge = -1.0 * coulombs;
    Spin neutronSpin = 0.5;
    Parity neutronParity = +1;

    THEN( "an exception is thrown" ) {

      // negative mass
      CHECK_THROWS( Particle( neutronID, negativeMass, neutronCharge,
                              neutronSpin, neutronParity ) );

      // negative charge
      CHECK_THROWS( Particle( neutronID, neutronMass, negativeCharge,
                              neutronSpin, neutronParity ) );
    } // THEN
  } // GIVEN
} // SCENARIO
