#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction/rmatrix/ParticleChannelData.hpp"

using namespace njoy::resonanceReconstruction;
using namespace dimwits::constant;

// convenience typedefs
using Particle = rmatrix::Particle;
using ParticlePair = rmatrix::ParticlePair;
using ParticleID = rmatrix::ParticleID;
using ParticlePairID = rmatrix::ParticlePairID;
using Neutron = rmatrix::Neutron;
using Photon = rmatrix::Photon;
using Fission = rmatrix::Fission;
using ChargedParticle = rmatrix::ChargedParticle;
template < typename Type > using Channel = rmatrix::Channel< Type >;
using ParticleChannel = rmatrix::ParticleChannel;
using ParticleChannelData = rmatrix::ParticleChannelData;
using ChannelID = rmatrix::ChannelID;
using ChannelQuantumNumbers = rmatrix::ChannelQuantumNumbers;
using ChannelRadii = rmatrix::ChannelRadii;
using BoundaryCondition = rmatrix::BoundaryCondition;

constexpr AtomicMass neutronMass = 1.008664 * daltons;
constexpr ElectricalCharge elementary = dimwits::constant::elementaryCharge;
constexpr double e = 1.6021766208e-19;

SCENARIO( "ParticleChannelData" ) {

  GIVEN( "valid data for a ParticleChannelData" ) {

    // particles
    Particle neutron( ParticleID( "n" ), 1.00866491582 * daltons,
                      0.0 * coulombs, 0.5, +1);
    Particle cl35( ParticleID( "Cl35" ), 34.968852694 * daltons,
                   17.0 * elementary, 1.5, +1);

    // particle pairs
    ParticlePair in( neutron, cl35 );

    // channel radii
    ChannelRadii radii( 4.822220e-1 * rootBarn,
                        3.667980e-1 * rootBarn );

    // channels
    Channel< Neutron > elastic( in, in, 0. * electronVolt,
                                { 0, 1.0, 1.0, +1 }, radii );

    // width and energy data
    std::vector< Energy > energies = { 1.2345e+1 * electronVolt,
                                       3500. * electronVolt };
    std::vector< ReducedWidth > widths = { 3.4e-2 * rootElectronVolt,
                                           6.5e-1 * rootElectronVolt };

    Energy energy = 1e-5 * electronVolt;

    THEN( "a ParticleChannelData instance can be constructed" ) {

      ParticleChannelData data( elastic,
                                std::move( energies ), std::move( widths ) );

      // verify channel data
      CHECK( "n,Cl35{0,1,1+}" == data.channelID() );
      CHECK( "n,Cl35->n,Cl35" == data.reactionID().symbol() );

      CHECK( true == data.isIncidentChannel() );
      CHECK( false == data.isEliminatedChannel() );

      auto incident = data.incidentParticlePair();
      CHECK( 1.00866491582 == Approx( incident.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident.particle().spin() ) );
      CHECK( +1 == incident.particle().parity() );
      CHECK( 34.968852694 == Approx( incident.residual().mass().value ) );
      CHECK( 17. * e == Approx( incident.residual().charge().value ) );
      CHECK( 1.5 == Approx( incident.residual().spin() ) );
      CHECK( +1 == incident.residual().parity() );
      CHECK( "n,Cl35" == incident.pairID().symbol() );

      auto pair = data.particlePair();
      CHECK( 1.00866491582 == Approx( pair.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair.particle().spin() ) );
      CHECK( +1 == pair.particle().parity() );
      CHECK( 34.968852694 == Approx( pair.residual().mass().value ) );
      CHECK( 17. * e == Approx( pair.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair.residual().spin() ) );
      CHECK( +1 == pair.residual().parity() );
      CHECK( "n,Cl35" == pair.pairID().symbol() );

      auto numbers = data.quantumNumbers();
      CHECK( 0 == numbers.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers.spin() );
      CHECK( 1.0 == numbers.totalAngularMomentum() );
      CHECK( +1 == numbers.parity() );

      auto radii = data.radii();
      CHECK( .4822220 == Approx( radii.penetrabilityRadius( energy ).value ) );
      CHECK( .4822220 == Approx( radii.shiftFactorRadius( energy ).value ) );
      CHECK( .3667980 == Approx( radii.phaseShiftRadius( energy ).value ) );

      CHECK( 0.0 == data.boundaryCondition() );

      CHECK( 0.0 == Approx( data.Q().value ) );

      // verify width data
      CHECK( 2 == data.energies().size() );
      CHECK( 12.345 == Approx( data.energies()[0].value ) );
      CHECK( 3500. == Approx( data.energies()[1].value ) );
      CHECK( 2 == data.widths().size() );
      CHECK( 3.4e-2 == Approx( data.widths()[0].value ) );
      CHECK( 0.65 == Approx( data.widths()[1].value ) );
    } // THEN
  } // GIVEN
} // SCENARIO
