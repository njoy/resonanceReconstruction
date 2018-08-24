#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using Particle = rmatrix::Particle;
using ParticlePair = rmatrix::ParticlePair;
using Neutron = rmatrix::Neutron;
using Fission = rmatrix::Fission;
using ChargedParticle = rmatrix::ChargedParticle;
using Photon = rmatrix::Photon;
template < typename Type > using Channel = rmatrix::Channel< Type >;
using Resonance = rmatrix::Resonance;
using ResonanceTable = rmatrix::ResonanceTable;
template < typename Option > using SpinGroup = rmatrix::SpinGroup< Option >;
using ShiftFactor = rmatrix::ShiftFactor;
using Constant = rmatrix::Constant;
using ParticleChannel = rmatrix::ParticleChannel;

constexpr AtomicMass neutronMass = 1.008664 * daltons;

SCENARIO( "SpinGroup" ) {

  GIVEN( "valid data for a SpinGroup" ) {

    // particles
    Particle photon( "g", 0.0 * daltons, 0.0 * coulombs, 1., +1);
    Particle neutron( "n", neutronMass, 0.0 * coulombs, 0.5, +1);
    Particle proton( "p", 1.007276 * daltons, 1.0 * coulombs, 0.5, +1);
    Particle pu240( "Pu240_e0", 2.379916e+2 * neutronMass,
                                94.0 * coulombs, 0.5, +1);
    Particle pu239( "Pu239_e0", 2.369986e+2 * neutronMass,
                                94.0 * coulombs, 0.5, +1);
    Particle np240( "Np240_e0", 2.379940e+2 * neutronMass,
                                93.0 * coulombs, 0.5, +1);
    Particle fission( "fission", 0.0 * daltons, 0.0 * coulombs, 0.0, +1);

    // particle pairs
    ParticlePair in( neutron, pu239, 0.0 * electronVolt );
    ParticlePair out1( photon, pu240, 0.0 * electronVolt );
    ParticlePair out2( neutron, fission, 0.0 * electronVolt, "fission" );
    ParticlePair out3( proton, np240, 0.0 * electronVolt );

    // channels
    Channel< Photon > capture( out1, { 0, 0.0, 1.0, +1 },
                               { 0.0 * rootBarn },
                               0.0 );
    Channel< Neutron > elastic( in, { 0, 0.5, 1.0, +1 },
                                { 9.410000e-1 * rootBarn },
                                0.0 );
    Channel< Fission > fission1( out2, "fission1", { 0, 0.0, 1.0, +1 },
                                 { 9.410000e-1 * rootBarn },
                                 0.0 );
    Channel< Fission > fission2( out2, "fission2", { 0, 0.0, 1.0, +1 },
                                 { 9.410000e-1 * rootBarn },
                                 0.0 );
    Channel< ChargedParticle > emission( out3, { 0, 0.5, 1.0, +1 },
                                         { 9.410000e-1 * rootBarn },
                                         0.0 );

    // single resonance table
    ResonanceTable single(
      { elastic.channelID(), fission1.channelID(),
        fission2.channelID(), emission.channelID() },
      { Resonance( 0.25 * electronVolt,
                   { 1.0 * rootElectronVolt, 2.0 * rootElectronVolt,
                     3.0 * rootElectronVolt, 4.0 * rootElectronVolt },
                   0.5 * rootElectronVolt ) } );

    THEN( "a SpinGroup can be constructed using the actual indices" ) {
      SpinGroup< ReichMoore, ShiftFactor >
          group( { 0 }, { elastic, fission1, fission2, emission },
                 std::move( single ) );

      REQUIRE( 1 == group.incidentChannels().size() );
      REQUIRE( "n,Pu239_e0{0,1/2,1+}" ==
                   std::visit( [] ( const auto& channel )
                                  { return channel.channelID(); },
                               group.incidentChannels().front() ) );

      REQUIRE( in.pairID() == group.incidentPair().pairID() );

      REQUIRE( 4 == group.channels().size() );

      // channel 1 - elastic
      auto channel1 =
          std::experimental::get< Channel< Neutron > >( group.channels()[0] );
      REQUIRE( "n,Pu239_e0{0,1/2,1+}" == channel1.channelID() );

      auto particlePair = channel1.particlePair();
      REQUIRE( 1.008664 == Approx( particlePair.particle().mass().value ) );
      REQUIRE( 0.0 == Approx( particlePair.particle().charge().value ) );
      REQUIRE( 0.5 == Approx( particlePair.particle().spin() ) );
      REQUIRE( +1 == particlePair.particle().parity() );
      REQUIRE( 239.0519559 == Approx( particlePair.residual().mass().value ) );
      REQUIRE( 94. == Approx( particlePair.residual().charge().value ) );
      REQUIRE( 0.5 == Approx( particlePair.residual().spin() ) );
      REQUIRE( +1 == particlePair.residual().parity() );
      REQUIRE( 0.0 == Approx( particlePair.Q().value ) );
      REQUIRE( "n,Pu239_e0" == particlePair.pairID() );

      auto numbers = channel1.quantumNumbers();
      REQUIRE( 0 == numbers.orbitalAngularMomentum() );
      REQUIRE( 0.5 == numbers.spin() );
      REQUIRE( 1.0 == numbers.totalAngularMomentum() );
      REQUIRE( +1 == numbers.parity() );

      auto radii = channel1.radii();
      REQUIRE( 0.941 ==
          Approx( radii.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      REQUIRE( 0.941 ==
          Approx( radii.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      REQUIRE( 0.941 ==
          Approx( radii.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      REQUIRE( 0.0 == channel1.boundaryCondition() );

      // channel 2 - fission1
      auto channel2 =
          std::experimental::get< Channel< Fission > >( group.channels()[1] );
      REQUIRE( "fission1{0,0,1+}" == channel2.channelID() );

      particlePair = channel2.particlePair();
      REQUIRE( 1.008664 == Approx( particlePair.particle().mass().value ) );
      REQUIRE( 0.0 == Approx( particlePair.particle().charge().value ) );
      REQUIRE( 0.5 == Approx( particlePair.particle().spin() ) );
      REQUIRE( +1 == particlePair.particle().parity() );
      REQUIRE( 0. == Approx( particlePair.residual().mass().value ) );
      REQUIRE( 0. == Approx( particlePair.residual().charge().value ) );
      REQUIRE( 0.0 == Approx( particlePair.residual().spin() ) );
      REQUIRE( +1 == particlePair.residual().parity() );
      REQUIRE( 0.0 == Approx( particlePair.Q().value ) );
      REQUIRE( "fission" == particlePair.pairID() );

      numbers = channel2.quantumNumbers();
      REQUIRE( 0 == numbers.orbitalAngularMomentum() );
      REQUIRE( 0.0 == numbers.spin() );
      REQUIRE( 1.0 == numbers.totalAngularMomentum() );
      REQUIRE( +1 == numbers.parity() );

      radii = channel2.radii();
      REQUIRE( 0.941 ==
          Approx( radii.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      REQUIRE( 0.941 ==
          Approx( radii.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      REQUIRE( 0.941 ==
          Approx( radii.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      REQUIRE( 0.0 == channel2.boundaryCondition() );

      // channel 3 - fission2
      auto channel3 =
          std::experimental::get< Channel< Fission > >( group.channels()[2] );
      REQUIRE( "fission2{0,0,1+}" == channel3.channelID() );

      particlePair = channel3.particlePair();
      REQUIRE( 1.008664 == Approx( particlePair.particle().mass().value ) );
      REQUIRE( 0.0 == Approx( particlePair.particle().charge().value ) );
      REQUIRE( 0.5 == Approx( particlePair.particle().spin() ) );
      REQUIRE( +1 == particlePair.particle().parity() );
      REQUIRE( 0. == Approx( particlePair.residual().mass().value ) );
      REQUIRE( 0. == Approx( particlePair.residual().charge().value ) );
      REQUIRE( 0.0 == Approx( particlePair.residual().spin() ) );
      REQUIRE( +1 == particlePair.residual().parity() );
      REQUIRE( 0.0 == Approx( particlePair.Q().value ) );
      REQUIRE( "fission" == particlePair.pairID() );

      numbers = channel3.quantumNumbers();
      REQUIRE( 0 == numbers.orbitalAngularMomentum() );
      REQUIRE( 0.0 == numbers.spin() );
      REQUIRE( 1.0 == numbers.totalAngularMomentum() );
      REQUIRE( +1 == numbers.parity() );

      radii = channel3.radii();
      REQUIRE( 0.941 ==
          Approx( radii.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      REQUIRE( 0.941 ==
          Approx( radii.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      REQUIRE( 0.941 ==
          Approx( radii.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      REQUIRE( 0.0 == channel3.boundaryCondition() );

      // channel 4 - emission
      auto channel4 =
          std::experimental::get< Channel< ChargedParticle > >
              ( group.channels()[3] );
      REQUIRE( "p,Np240_e0{0,1/2,1+}" == channel4.channelID() );

      particlePair = channel4.particlePair();
      REQUIRE( 1.007276 == Approx( particlePair.particle().mass().value ) );
      REQUIRE( 1.0 == Approx( particlePair.particle().charge().value ) );
      REQUIRE( 0.5 == Approx( particlePair.particle().spin() ) );
      REQUIRE( +1 == particlePair.particle().parity() );
      REQUIRE( 240.055980 == Approx( particlePair.residual().mass().value ) );
      REQUIRE( 93. == Approx( particlePair.residual().charge().value ) );
      REQUIRE( 0.5 == Approx( particlePair.residual().spin() ) );
      REQUIRE( +1 == particlePair.residual().parity() );
      REQUIRE( 0.0 == Approx( particlePair.Q().value ) );
      REQUIRE( "p,Np240_e0" == particlePair.pairID() );

      numbers = channel4.quantumNumbers();
      REQUIRE( 0 == numbers.orbitalAngularMomentum() );
      REQUIRE( 0.5 == numbers.spin() );
      REQUIRE( 1.0 == numbers.totalAngularMomentum() );
      REQUIRE( +1 == numbers.parity() );

      radii = channel4.radii();
      REQUIRE( 0.941 ==
          Approx( radii.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      REQUIRE( 0.941 ==
          Approx( radii.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      REQUIRE( 0.941 ==
          Approx( radii.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      REQUIRE( 0.0 == channel4.boundaryCondition() );

      // reactions in the spin group
      auto reactions = group.reactions();
      REQUIRE( 5 == reactions.size() );
      REQUIRE( "n,Pu239_e0->n,Pu239_e0" == reactions[0] );
      REQUIRE( "n,Pu239_e0->fission" == reactions[1] );
      REQUIRE( "n,Pu239_e0->fission" == reactions[2] );
      REQUIRE( "n,Pu239_e0->p,Np240_e0" == reactions[3] );
      REQUIRE( "n,Pu239_e0->capture" == reactions[4] );

      // resonance data
      auto table = group.resonanceTable();
      REQUIRE( 4 == table.numberChannels() );
      REQUIRE( 4 == table.channels().size() );
      REQUIRE( "n,Pu239_e0{0,1/2,1+}" == table.channels()[0] );
      REQUIRE( "fission1{0,0,1+}" == table.channels()[1] );
      REQUIRE( "fission2{0,0,1+}" == table.channels()[2] );
      REQUIRE( "p,Np240_e0{0,1/2,1+}" == table.channels()[3] );
      REQUIRE( 1 == table.numberResonances() );
      REQUIRE( 1 == table.resonances().size() );
      REQUIRE( 1 == table.energies().size() );
      REQUIRE( 0.25 == Approx( table.energies()[0].value ) );
      auto resonance = table.resonances()[0];
      REQUIRE( 0.25 == Approx( resonance.energy().value ) );
      REQUIRE( 0.5 == Approx( resonance.eliminatedWidth().value ) );
      REQUIRE( 4 == resonance.widths().size() );
      REQUIRE( 1.0 == Approx( resonance.widths()[0].value ) );
      REQUIRE( 2.0 == Approx( resonance.widths()[1].value ) );
      REQUIRE( 3.0 == Approx( resonance.widths()[2].value ) );
      REQUIRE( 4.0 == Approx( resonance.widths()[3].value ) );
    }

    THEN( "a SpinGroup can be constructed using the incident particle pair" ) {
      SpinGroup< ReichMoore, ShiftFactor >
          group( in, { elastic, fission1, fission2, emission },
                 std::move( single ) );

      REQUIRE( 1 == group.incidentChannels().size() );
      REQUIRE( "n,Pu239_e0{0,1/2,1+}" ==
                   std::visit( [] ( const auto& channel )
                                  { return channel.channelID(); },
                               group.incidentChannels().front() ) );

      REQUIRE( in.pairID() == group.incidentPair().pairID() );

      REQUIRE( 4 == group.channels().size() );

      // channel 1 - elastic
      auto channel1 =
          std::experimental::get< Channel< Neutron > >( group.channels()[0] );
      REQUIRE( "n,Pu239_e0{0,1/2,1+}" == channel1.channelID() );

      auto particlePair = channel1.particlePair();
      REQUIRE( 1.008664 == Approx( particlePair.particle().mass().value ) );
      REQUIRE( 0.0 == Approx( particlePair.particle().charge().value ) );
      REQUIRE( 0.5 == Approx( particlePair.particle().spin() ) );
      REQUIRE( +1 == particlePair.particle().parity() );
      REQUIRE( 239.0519559 == Approx( particlePair.residual().mass().value ) );
      REQUIRE( 94. == Approx( particlePair.residual().charge().value ) );
      REQUIRE( 0.5 == Approx( particlePair.residual().spin() ) );
      REQUIRE( +1 == particlePair.residual().parity() );
      REQUIRE( 0.0 == Approx( particlePair.Q().value ) );
      REQUIRE( "n,Pu239_e0" == particlePair.pairID() );

      auto numbers = channel1.quantumNumbers();
      REQUIRE( 0 == numbers.orbitalAngularMomentum() );
      REQUIRE( 0.5 == numbers.spin() );
      REQUIRE( 1.0 == numbers.totalAngularMomentum() );
      REQUIRE( +1 == numbers.parity() );

      auto radii = channel1.radii();
      REQUIRE( 0.941 ==
          Approx( radii.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      REQUIRE( 0.941 ==
          Approx( radii.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      REQUIRE( 0.941 ==
          Approx( radii.phaseShiftRadius( 1e-5 * electronVolt ).value ) );
      REQUIRE( 0.0 == channel1.boundaryCondition() );

      // channel 2 - fission1
      auto channel2 =
          std::experimental::get< Channel< Fission > >( group.channels()[1] );
      REQUIRE( "fission1{0,0,1+}" == channel2.channelID() );

      particlePair = channel2.particlePair();
      REQUIRE( 1.008664 == Approx( particlePair.particle().mass().value ) );
      REQUIRE( 0.0 == Approx( particlePair.particle().charge().value ) );
      REQUIRE( 0.5 == Approx( particlePair.particle().spin() ) );
      REQUIRE( +1 == particlePair.particle().parity() );
      REQUIRE( 0. == Approx( particlePair.residual().mass().value ) );
      REQUIRE( 0. == Approx( particlePair.residual().charge().value ) );
      REQUIRE( 0.0 == Approx( particlePair.residual().spin() ) );
      REQUIRE( +1 == particlePair.residual().parity() );
      REQUIRE( 0.0 == Approx( particlePair.Q().value ) );
      REQUIRE( "fission" == particlePair.pairID() );

      numbers = channel2.quantumNumbers();
      REQUIRE( 0 == numbers.orbitalAngularMomentum() );
      REQUIRE( 0.0 == numbers.spin() );
      REQUIRE( 1.0 == numbers.totalAngularMomentum() );
      REQUIRE( +1 == numbers.parity() );

      radii = channel2.radii();
      REQUIRE( 0.941 ==
          Approx( radii.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      REQUIRE( 0.941 ==
          Approx( radii.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      REQUIRE( 0.941 ==
          Approx( radii.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      REQUIRE( 0.0 == channel2.boundaryCondition() );

      // channel 3 - fission2
      auto channel3 =
          std::experimental::get< Channel< Fission > >( group.channels()[2] );
      REQUIRE( "fission2{0,0,1+}" == channel3.channelID() );

      particlePair = channel3.particlePair();
      REQUIRE( 1.008664 == Approx( particlePair.particle().mass().value ) );
      REQUIRE( 0.0 == Approx( particlePair.particle().charge().value ) );
      REQUIRE( 0.5 == Approx( particlePair.particle().spin() ) );
      REQUIRE( +1 == particlePair.particle().parity() );
      REQUIRE( 0. == Approx( particlePair.residual().mass().value ) );
      REQUIRE( 0. == Approx( particlePair.residual().charge().value ) );
      REQUIRE( 0.0 == Approx( particlePair.residual().spin() ) );
      REQUIRE( +1 == particlePair.residual().parity() );
      REQUIRE( 0.0 == Approx( particlePair.Q().value ) );
      REQUIRE( "fission" == particlePair.pairID() );

      numbers = channel3.quantumNumbers();
      REQUIRE( 0 == numbers.orbitalAngularMomentum() );
      REQUIRE( 0.0 == numbers.spin() );
      REQUIRE( 1.0 == numbers.totalAngularMomentum() );
      REQUIRE( +1 == numbers.parity() );

      radii = channel3.radii();
      REQUIRE( 0.941 ==
          Approx( radii.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      REQUIRE( 0.941 ==
          Approx( radii.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      REQUIRE( 0.941 ==
          Approx( radii.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      REQUIRE( 0.0 == channel3.boundaryCondition() );

      // channel 4 - emission
      auto channel4 =
          std::experimental::get< Channel< ChargedParticle > >
              ( group.channels()[3] );
      REQUIRE( "p,Np240_e0{0,1/2,1+}" == channel4.channelID() );

      particlePair = channel4.particlePair();
      REQUIRE( 1.007276 == Approx( particlePair.particle().mass().value ) );
      REQUIRE( 1.0 == Approx( particlePair.particle().charge().value ) );
      REQUIRE( 0.5 == Approx( particlePair.particle().spin() ) );
      REQUIRE( +1 == particlePair.particle().parity() );
      REQUIRE( 240.055980 == Approx( particlePair.residual().mass().value ) );
      REQUIRE( 93. == Approx( particlePair.residual().charge().value ) );
      REQUIRE( 0.5 == Approx( particlePair.residual().spin() ) );
      REQUIRE( +1 == particlePair.residual().parity() );
      REQUIRE( 0.0 == Approx( particlePair.Q().value ) );
      REQUIRE( "p,Np240_e0" == particlePair.pairID() );

      numbers = channel4.quantumNumbers();
      REQUIRE( 0 == numbers.orbitalAngularMomentum() );
      REQUIRE( 0.5 == numbers.spin() );
      REQUIRE( 1.0 == numbers.totalAngularMomentum() );
      REQUIRE( +1 == numbers.parity() );

      radii = channel4.radii();
      REQUIRE( 0.941 ==
          Approx( radii.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      REQUIRE( 0.941 ==
          Approx( radii.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      REQUIRE( 0.941 ==
          Approx( radii.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      REQUIRE( 0.0 == channel4.boundaryCondition() );

      // reactions in the spin group
      auto reactions = group.reactions();
      REQUIRE( 5 == reactions.size() );
      REQUIRE( "n,Pu239_e0->n,Pu239_e0" == reactions[0] );
      REQUIRE( "n,Pu239_e0->fission" == reactions[1] );
      REQUIRE( "n,Pu239_e0->fission" == reactions[2] );
      REQUIRE( "n,Pu239_e0->p,Np240_e0" == reactions[3] );
      REQUIRE( "n,Pu239_e0->capture" == reactions[4] );

      // resonance data
      auto table = group.resonanceTable();
      REQUIRE( 4 == table.numberChannels() );
      REQUIRE( 4 == table.channels().size() );
      REQUIRE( "n,Pu239_e0{0,1/2,1+}" == table.channels()[0] );
      REQUIRE( "fission1{0,0,1+}" == table.channels()[1] );
      REQUIRE( "fission2{0,0,1+}" == table.channels()[2] );
      REQUIRE( "p,Np240_e0{0,1/2,1+}" == table.channels()[3] );
      REQUIRE( 1 == table.numberResonances() );
      REQUIRE( 1 == table.resonances().size() );
      REQUIRE( 1 == table.energies().size() );
      REQUIRE( 0.25 == Approx( table.energies()[0].value ) );
      auto resonance = table.resonances()[0];
      REQUIRE( 0.25 == Approx( resonance.energy().value ) );
      REQUIRE( 0.5 == Approx( resonance.eliminatedWidth().value ) );
      REQUIRE( 4 == resonance.widths().size() );
      REQUIRE( 1.0 == Approx( resonance.widths()[0].value ) );
      REQUIRE( 2.0 == Approx( resonance.widths()[1].value ) );
      REQUIRE( 3.0 == Approx( resonance.widths()[2].value ) );
      REQUIRE( 4.0 == Approx( resonance.widths()[3].value ) );
    }
  } // GIVEN

  GIVEN( "data for a SpinGroup with incoherent J,pi channels" ) {

    // particles
    Particle photon( "g", 0.0 * daltons, 0.0 * coulombs, 1., +1);
    Particle neutron( "n", neutronMass, 0.0 * coulombs, 0.5, +1);
    Particle proton( "p", 1.007276 * daltons, 1.0 * coulombs, 0.5, +1);
    Particle pu240( "Pu240_e0", 2.379916e+2 * neutronMass,
                                94.0 * coulombs, 0.5, +1);
    Particle pu239( "Pu239_e0", 2.369986e+2 * neutronMass,
                                94.0 * coulombs, 0.5, +1);
    Particle np240( "Np240_e0", 2.379940e+2 * neutronMass,
                                93.0 * coulombs, 0.5, +1);
    Particle fission( "fission", 0.0 * daltons, 0.0 * coulombs, 0.0, +1);

    // particle pairs
    ParticlePair in( neutron, pu239, 0.0 * electronVolt );
    ParticlePair out1( photon, pu240, 0.0 * electronVolt );
    ParticlePair out2( neutron, fission, 0.0 * electronVolt, "fission" );
    ParticlePair out3( proton, np240, 0.0 * electronVolt );

    // channels (the second channel has an incorrect J,pi with respect to the
    // other channels)
    Channel< Photon > capture( out1, { 0, 0.0, 1.0, +1 },
                               { 0.0 * rootBarn },
                               0.0 );
    Channel< Neutron > incorrect( in, { 0, 0.5, 0.0, +1 },
                                  { 9.410000e-1 * rootBarn },
                                  0.0 );
    Channel< Fission > fission1( out2, "fission1", { 0, 0.0, 1.0, +1 },
                                 { 9.410000e-1 * rootBarn },
                                 0.0 );
    Channel< Fission > fission2( out2, "fission2", { 0, 0.0, 1.0, +1 },
                                 { 9.410000e-1 * rootBarn },
                                 0.0 );
    Channel< ChargedParticle > emission( out3, { 0, 0.5, 1.0, +1 },
                                         { 9.410000e-1 * rootBarn },
                                         0.0 );

    // single resonance table
    ResonanceTable single(
      { incorrect.channelID(), fission1.channelID(),
        fission2.channelID(), emission.channelID() },
      { Resonance( 0.25 * electronVolt,
                   { 1.0 * rootElectronVolt, 2.0 * rootElectronVolt,
                     3.0 * rootElectronVolt, 4.0 * rootElectronVolt },
                   0.5 * rootElectronVolt ) } );

    THEN( "an exception is thrown at construction" ) {

      REQUIRE_THROWS( SpinGroup< ReichMoore, ShiftFactor >
                          ( in, { incorrect, fission1, fission2, emission },
                            std::move( single ) ) );
    }
  } // GIVEN
} // SCENARIO


