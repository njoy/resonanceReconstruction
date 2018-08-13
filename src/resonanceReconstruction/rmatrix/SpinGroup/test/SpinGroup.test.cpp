#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using Particle = rmatrix::Particle;
using ParticlePair = rmatrix::ParticlePair;
using Neutron = rmatrix::Neutron;
using Photon = rmatrix::Photon;
template < typename Type > using Channel = rmatrix::Channel< Type >;
using Resonance = rmatrix::Resonance;
using ResonanceTable = rmatrix::ResonanceTable;
template < typename Option > using SpinGroup = rmatrix::SpinGroup< Option >;
using Sammy = rmatrix::Sammy;
using Constant = rmatrix::Constant;

constexpr AtomicMass neutronMass = 1.008664 * daltons;

SCENARIO( "SpinGroup" ) {

  GIVEN( "valid data for a SpinGroup" ) {

    // test based on Fe54 ENDF/B-VIII.0 LRF7 resonance evaluation
    // data given in Gamma = 2 gamma^2 P(Er) so conversion is required

    Particle neutron( "n", neutronMass, 0.0 * coulombs, 0.5, +1);
    Particle fe54( "Fe54_e0", 54. * daltons, 26.0 * coulombs, 0.0, +1);
    ParticlePair pair( neutron, fe54, 0.0 * electronVolt );
    Channel< Neutron > elastic( "2", pair, { 0, 0.5, 0.5, +1 },
                                { 5.437300e-1 * rootBarn,
                                  5.437300e-1 * rootBarn },
                                0.0 );

    // single resonance table
    ResonanceTable single(
      { "2" },
      { Resonance( 1.0 * electronVolt, { 2.0 * rootElectronVolt },
                   0.5 * rootElectronVolt ) } );

    THEN( "a SpinGroup can be constructed" ) {
      SpinGroup< Sammy > group( { elastic }, { 0 }, std::move( single ) );

      REQUIRE( 1 == group.channels().size() );

      auto channel = std::experimental::get< Channel< Neutron > >( group.channels()[0] );
      REQUIRE( "2" == channel.channelID() );
      REQUIRE( 1.008664 == Approx( channel.particlePair().particle().mass().value ) );
      REQUIRE( 0.0 == Approx( channel.particlePair().particle().charge().value ) );
      REQUIRE( 0.5 == Approx( channel.particlePair().particle().spin() ) );
      REQUIRE( +1 == channel.particlePair().particle().parity() );
      REQUIRE( 54. == Approx( channel.particlePair().residual().mass().value ) );
      REQUIRE( 26. == Approx( channel.particlePair().residual().charge().value ) );
      REQUIRE( 0.0 == Approx( channel.particlePair().residual().spin() ) );
      REQUIRE( +1 == channel.particlePair().residual().parity() );
      REQUIRE( 0.0 == Approx( channel.particlePair().Q().value ) );
      REQUIRE( "n,Fe54_e0" == channel.particlePair().pairID() );
      REQUIRE( 0 == channel.quantumNumbers().orbitalAngularMomentum() );
      REQUIRE( 0.5 == channel.quantumNumbers().spin() );
      REQUIRE( 0.5 == channel.quantumNumbers().totalAngularMomentum() );
      REQUIRE( +1 == channel.quantumNumbers().parity() );
      REQUIRE( 5.437300e-1 ==
          Approx( channel.radii().penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      REQUIRE( 5.437300e-1 ==
          Approx( channel.radii().shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      REQUIRE( 5.437300e-1 ==
          Approx( channel.radii().phaseShiftRadius( 1e-5 * electronVolt ).value ) );
      REQUIRE( 0.0 == channel.boundaryCondition() );

      REQUIRE( 2 == group.reactions().size() );
      REQUIRE( "n,Fe54_e0->n,Fe54_e0" == group.reactions()[0] );
      REQUIRE( "n,Fe54_e0->capture" == group.reactions()[1] );

      auto table = group.resonanceTable();
      REQUIRE( 1 == table.numberChannels() );
      REQUIRE( 1 == table.channels().size() );
      REQUIRE( "2" == table.channels()[0] );
      REQUIRE( 1 == table.numberResonances() );
      REQUIRE( 1 == table.resonances().size() );
      REQUIRE( 1 == table.energies().size() );
      REQUIRE( 1.0 == Approx( table.energies()[0].value ) );
      auto resonance = table.resonances()[0];
      REQUIRE( 1.0 == Approx( resonance.energy().value ) );
      REQUIRE( 0.5 == Approx( resonance.eliminatedWidth().value ) );
      REQUIRE( 1 == resonance.widths().size() );
      REQUIRE( 2.0 == Approx( resonance.widths()[0].value ) );
    }
  } // GIVEN
} // SCENARIO


