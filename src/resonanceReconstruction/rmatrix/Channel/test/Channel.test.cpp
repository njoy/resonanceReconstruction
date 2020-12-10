#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction/rmatrix/Channel.hpp"

using namespace njoy::resonanceReconstruction;

// convenience typedefs
using Particle = rmatrix::Particle;
using ParticleID = rmatrix::ParticleID;
using ParticlePair = rmatrix::ParticlePair;
using ParticlePairID = rmatrix::ParticlePairID;
using Neutron = rmatrix::Neutron;
using Photon = rmatrix::Photon;
using ChargedParticle = rmatrix::ChargedParticle;
template < typename Type > using Channel = rmatrix::Channel< Type >;
using ChannelID = rmatrix::ChannelID;
using ChannelQuantumNumbers = rmatrix::ChannelQuantumNumbers;
using ChannelRadii = rmatrix::ChannelRadii;
using BoundaryCondition = rmatrix::BoundaryCondition;

constexpr ElectricalCharge elementary = dimwits::constant::elementaryCharge;
constexpr double e = 1.6021766208e-19;

SCENARIO( "Channel" ) {

  GIVEN( "valid data for a Channel" ) {

    // particles
    Particle photon( ParticleID( "g" ), 0.0 * daltons, 0.0 * coulombs, 1., +1);
    Particle neutron( ParticleID( "n" ), 1.00866491582 * daltons,
                      0.0 * coulombs, 0.5, +1);
    Particle proton( ParticleID( "p" ), 1.00727647 * daltons, elementary, 0.5, +1);
    Particle cl36( ParticleID( "Cl36_e0" ), 35.968306822 * daltons,
                   17.0 * elementary, 0., +1);
    Particle cl35( ParticleID( "Cl35_e0" ), 34.968852694 * daltons,
                   17.0 * elementary, 1.5, +1);
    Particle cl35_e1( ParticleID( "Cl35_e1" ), 34.968852694 * daltons,
                      17.0 * elementary, 1.5, +1);
    Particle s36( ParticleID( "S36_e0" ), 35.967080699 * daltons,
                  16.0 * elementary, 1.5, +1);

    // particle pairs
    ParticlePair elasticPair( neutron, cl35 );
    ParticlePair inelasticPair( neutron, cl35_e1 );
    ParticlePair capturePair( photon, cl36 );
    ParticlePair protonEmissionPair( proton, s36 );

    // Q values
    QValue elasticQ = 0.0 * electronVolt;
    QValue inelasticQ = -1.219440e+6 * electronVolt;
    QValue captureQ = 0.0 * electronVolt;
    QValue protonEmissionQ = 6.152200e+5 * electronVolt;

    // quantum numbers
    ChannelQuantumNumbers elasticNumbers( 0, 1.0, 1.0, +1 );
    ChannelQuantumNumbers inelasticNumbers( 0, 1.0, 1.0, +1 );
    ChannelQuantumNumbers captureNumbers( 0, 0.0, 1.0, +1 );
    ChannelQuantumNumbers protonEmissionNumbers( 0, 1.0, 1.0, +1 );

    // channel radii
    ChannelRadii elasticRadii( 4.822220e-1 * rootBarn,
                               3.667980e-1 * rootBarn );
    ChannelRadii inelasticRadii( 4.822220e-1 * rootBarn,
                               3.667980e-1 * rootBarn );
    ChannelRadii captureRadii( 0.0 * rootBarn );
    ChannelRadii protonEmissionRadii( 4.822220e-1 * rootBarn,
                                      3.667980e-1 * rootBarn );

    // custom channel identifier
    ChannelID id = "1";

    THEN( "Channels can be constructed with an automatic ChannelID" ) {

      Energy energy = 1e-5 * electronVolt;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // photon channel, not an incident channel, no threshold
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      Channel< Photon > capture( elasticPair, capturePair, captureQ,
                                 captureNumbers, captureRadii );

      CHECK( "photon,Cl36{0,0,1+}" == capture.channelID() );
      CHECK( "n,Cl35->photon,Cl36" == capture.reactionID().symbol() );

      CHECK( false == capture.isIncidentChannel() );

      auto incident = capture.incidentParticlePair();
      CHECK( 1.00866491582 == Approx( incident.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident.particle().spin() ) );
      CHECK( +1 == incident.particle().parity() );
      CHECK( 34.968852694 == Approx( incident.residual().mass().value ) );
      CHECK( 17. * e == Approx( incident.residual().charge().value ) );
      CHECK( 1.5 == Approx( incident.residual().spin() ) );
      CHECK( +1 == incident.residual().parity() );
      CHECK( "n,Cl35" == incident.pairID().symbol() );

      auto pair = capture.particlePair();
      CHECK( 0.0 == Approx( pair.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair.particle().charge().value ) );
      CHECK( 1.0 == Approx( pair.particle().spin() ) );
      CHECK( +1 == pair.particle().parity() );
      CHECK( 35.968306822 == Approx( pair.residual().mass().value ) );
      CHECK( 17. * e == Approx( pair.residual().charge().value ) );
      CHECK( 0.0 == Approx( pair.residual().spin() ) );
      CHECK( +1 == pair.residual().parity() );
      CHECK( "photon,Cl36" == pair.pairID().symbol() );

      auto numbers = capture.quantumNumbers();
      CHECK( 0 == numbers.orbitalAngularMomentum() );
      CHECK( 0.0 == numbers.spin() );
      CHECK( 1.0 == numbers.totalAngularMomentum() );
      CHECK( +1 == numbers.parity() );

      auto radii = capture.radii();
      CHECK( 0.0 == Approx( radii.penetrabilityRadius( energy ).value ) );
      CHECK( 0.0 == Approx( radii.shiftFactorRadius( energy ).value ) );
      CHECK( 0.0 == Approx( radii.phaseShiftRadius( energy ).value ) );

      CHECK( 0.0 == capture.boundaryCondition() );

      CHECK( 0.0 == Approx( capture.Q().value ) );

      CHECK( 1.0 == Approx( capture.statisticalSpinFactor() ) );

      CHECK( false == Approx( capture.belowThreshold( energy ) ) );

      CHECK( 0.0 == Approx( capture.sommerfeldParameter( energy ) ) );
      CHECK( 0.0 == Approx( capture.waveNumber( energy ).value ) );
      CHECK( 1.0 == Approx( capture.penetrability( energy ) ) );
      CHECK( 0.0 == Approx( capture.shiftFactor( energy ) ) );
      CHECK( 0.0 == Approx( capture.phaseShift( energy ) ) );
      CHECK( 0.0 == Approx( capture.coulombPhaseShift( energy ) ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // neutron channel, incident channel, no threshold
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      Channel< Neutron > elastic( elasticPair, elasticPair, elasticQ,
                                  elasticNumbers, elasticRadii );

      CHECK( "n,Cl35{0,1,1+}" == elastic.channelID() );
      CHECK( "n,Cl35->n,Cl35" == elastic.reactionID().symbol() );

      CHECK( true == elastic.isIncidentChannel() );

      incident = elastic.incidentParticlePair();
      CHECK( 1.00866491582 == Approx( incident.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident.particle().spin() ) );
      CHECK( +1 == incident.particle().parity() );
      CHECK( 34.968852694 == Approx( incident.residual().mass().value ) );
      CHECK( 17. * e == Approx( incident.residual().charge().value ) );
      CHECK( 1.5 == Approx( incident.residual().spin() ) );
      CHECK( +1 == incident.residual().parity() );
      CHECK( "n,Cl35" == incident.pairID().symbol() );

      pair = elastic.particlePair();
      CHECK( 1.00866491582 == Approx( pair.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair.particle().spin() ) );
      CHECK( +1 == pair.particle().parity() );
      CHECK( 34.968852694 == Approx( pair.residual().mass().value ) );
      CHECK( 17. * e == Approx( pair.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair.residual().spin() ) );
      CHECK( +1 == pair.residual().parity() );
      CHECK( "n,Cl35" == pair.pairID().symbol() );

      numbers = elastic.quantumNumbers();
      CHECK( 0 == numbers.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers.spin() );
      CHECK( 1.0 == numbers.totalAngularMomentum() );
      CHECK( +1 == numbers.parity() );

      radii = elastic.radii();
      CHECK( .4822220 == Approx( radii.penetrabilityRadius( energy ).value ) );
      CHECK( .4822220 == Approx( radii.shiftFactorRadius( energy ).value ) );
      CHECK( .3667980 == Approx( radii.phaseShiftRadius( energy ).value ) );

      CHECK( 0.0 == elastic.boundaryCondition() );

      CHECK( 0.0 == Approx( elastic.Q().value ) );

      CHECK( 0.375 == Approx( elastic.statisticalSpinFactor() ) );

      CHECK( false == Approx( elastic.belowThreshold( energy ) ) );

      CHECK( 0.0 == Approx( elastic.sommerfeldParameter( energy ) ) );
      CHECK( 6.75215238E-06 == Approx( elastic.waveNumber( energy ).value ) );
      CHECK( 3.25603642E-06 == Approx( elastic.penetrability( energy ) ) );
      CHECK( 0.0 == Approx( elastic.shiftFactor( energy ) ) );
      CHECK( 2.47667599E-06 == Approx( elastic.phaseShift( energy ) ) );
      CHECK( 0.0 == Approx( elastic.coulombPhaseShift( energy ) ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // neutron channel, not incident channel, threshold
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      Channel< Neutron > inelastic( elasticPair, inelasticPair, inelasticQ,
                                    inelasticNumbers, inelasticRadii );

      CHECK( "n,Cl35_e1{0,1,1+}" == inelastic.channelID() );
      CHECK( "n,Cl35->n,Cl35_e1" == inelastic.reactionID().symbol() );

      CHECK( false == inelastic.isIncidentChannel() );

      incident = inelastic.incidentParticlePair();
      CHECK( 1.00866491582 == Approx( incident.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident.particle().spin() ) );
      CHECK( +1 == incident.particle().parity() );
      CHECK( 34.968852694 == Approx( incident.residual().mass().value ) );
      CHECK( 17. * e == Approx( incident.residual().charge().value ) );
      CHECK( 1.5 == Approx( incident.residual().spin() ) );
      CHECK( +1 == incident.residual().parity() );
      CHECK( "n,Cl35" == incident.pairID().symbol() );

      pair = inelastic.particlePair();
      CHECK( 1.00866491582 == Approx( pair.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair.particle().spin() ) );
      CHECK( +1 == pair.particle().parity() );
      CHECK( 34.968852694 == Approx( pair.residual().mass().value ) );
      CHECK( 17. * e == Approx( pair.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair.residual().spin() ) );
      CHECK( +1 == pair.residual().parity() );
      CHECK( "n,Cl35_e1" == pair.pairID().symbol() );

      numbers = inelastic.quantumNumbers();
      CHECK( 0 == numbers.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers.spin() );
      CHECK( 1.0 == numbers.totalAngularMomentum() );
      CHECK( +1 == numbers.parity() );

      radii = inelastic.radii();
      CHECK( .4822220 == Approx( radii.penetrabilityRadius( energy ).value ) );
      CHECK( .4822220 == Approx( radii.shiftFactorRadius( energy ).value ) );
      CHECK( .3667980 == Approx( radii.phaseShiftRadius( energy ).value ) );

      CHECK( 0.0 == inelastic.boundaryCondition() );

      CHECK( -1.219440e+6 == Approx( inelastic.Q().value ) );

      CHECK( 0.375 == Approx( inelastic.statisticalSpinFactor() ) );

      CHECK( true == Approx( inelastic.belowThreshold( energy ) ) );
      CHECK( false
             == Approx( inelastic.belowThreshold( 15e+6 * electronVolts ) ) );

      CHECK( 0.0 == Approx( inelastic.sommerfeldParameter( energy ) ) );
      CHECK( 2.39164854E+00 == Approx( inelastic.waveNumber( energy ).value ) );
      CHECK( 1.15330554E+00 == Approx( inelastic.penetrability( energy ) ) );
      CHECK( 0.0 == Approx( inelastic.shiftFactor( energy ) ) );
      CHECK( 8.77251899E-01 == Approx( inelastic.phaseShift( energy ) ) );
      CHECK( 0.0 == Approx( inelastic.coulombPhaseShift( energy ) ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // proton channel, not an incident channel, no threshold
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      Channel< ChargedParticle > protonEmission( elasticPair,
                                                 protonEmissionPair,
                                                 protonEmissionQ,
                                                 protonEmissionNumbers,
                                                 protonEmissionRadii );

      CHECK( "p,S36{0,1,1+}" == protonEmission.channelID() );
      CHECK( "n,Cl35->p,S36" == protonEmission.reactionID().symbol() );

      CHECK( false == protonEmission.isIncidentChannel() );

      incident = protonEmission.incidentParticlePair();
      CHECK( 1.00866491582 == Approx( incident.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident.particle().spin() ) );
      CHECK( +1 == incident.particle().parity() );
      CHECK( 34.968852694 == Approx( incident.residual().mass().value ) );
      CHECK( 17. * e == Approx( incident.residual().charge().value ) );
      CHECK( 1.5 == Approx( incident.residual().spin() ) );
      CHECK( +1 == incident.residual().parity() );
      CHECK( "n,Cl35" == incident.pairID().symbol() );

      pair = protonEmission.particlePair();
      CHECK( 1.00727647 == Approx( pair.particle().mass().value ) );
      CHECK( e == Approx( pair.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair.particle().spin() ) );
      CHECK( +1 == pair.particle().parity() );
      CHECK( 35.967080699 == Approx( pair.residual().mass().value ) );
      CHECK( 16. * e == Approx( pair.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair.residual().spin() ) );
      CHECK( +1 == pair.residual().parity() );
      CHECK( "p,S36" == pair.pairID().symbol() );

      numbers = protonEmission.quantumNumbers();
      CHECK( 0 == numbers.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers.spin() );
      CHECK( 1.0 == numbers.totalAngularMomentum() );
      CHECK( +1 == numbers.parity() );

      radii = protonEmission.radii();
      CHECK( .4822220 == Approx( radii.penetrabilityRadius( energy ).value ) );
      CHECK( .4822220 == Approx( radii.shiftFactorRadius( energy ).value ) );
      CHECK( .3667980 == Approx( radii.phaseShiftRadius( energy ).value ) );

      CHECK( 0.0 == protonEmission.boundaryCondition() );

      CHECK( 6.152200e+5 == Approx( protonEmission.Q().value ) );

      CHECK( 0.375 == Approx( protonEmission.statisticalSpinFactor() ) );

      CHECK( false == Approx( protonEmission.belowThreshold( energy ) ) );

      CHECK( 3.17996084E+00 == Approx( protonEmission.sommerfeldParameter( energy ) ) );
      CHECK( 1.69828445E+00 == Approx( protonEmission.waveNumber( energy ).value ) );
      CHECK( 0.000027793 == Approx( protonEmission.penetrability( energy ) ) );
      CHECK( -1.8734549658 == Approx( protonEmission.shiftFactor( energy ) ) );
      CHECK( 0.0000020468 == Approx( protonEmission.phaseShift( energy ) ) );
      CHECK( 0.0 == Approx( protonEmission.coulombPhaseShift( energy ) ) );

      // remark: P, S and phi for proton channel copied from code output
      // other values calculated by hand
    } // THEN

    THEN( "a Channel can be constructed with an alternate ParticleID" ) {

      Energy energy = 1e-5 * electronVolt;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // photon channel, not an incident channel
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      Channel< Photon > capture( elasticPair, capturePair, id, captureQ,
                                 captureNumbers, captureRadii );

      CHECK( "1{0,0,1+}" == capture.channelID() );
      CHECK( "n,Cl35->photon,Cl36" == capture.reactionID().symbol() );

      CHECK( false == capture.isIncidentChannel() );

      auto incident = capture.incidentParticlePair();
      CHECK( 1.00866491582 == Approx( incident.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident.particle().spin() ) );
      CHECK( +1 == incident.particle().parity() );
      CHECK( 34.968852694 == Approx( incident.residual().mass().value ) );
      CHECK( 17. * e == Approx( incident.residual().charge().value ) );
      CHECK( 1.5 == Approx( incident.residual().spin() ) );
      CHECK( +1 == incident.residual().parity() );
      CHECK( "n,Cl35" == incident.pairID().symbol() );

      auto pair = capture.particlePair();
      CHECK( 0.0 == Approx( pair.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair.particle().charge().value ) );
      CHECK( 1.0 == Approx( pair.particle().spin() ) );
      CHECK( +1 == pair.particle().parity() );
      CHECK( 35.968306822 == Approx( pair.residual().mass().value ) );
      CHECK( 17. * e == Approx( pair.residual().charge().value ) );
      CHECK( 0.0 == Approx( pair.residual().spin() ) );
      CHECK( +1 == pair.residual().parity() );
      CHECK( "photon,Cl36" == pair.pairID().symbol() );

      auto numbers = capture.quantumNumbers();
      CHECK( 0 == numbers.orbitalAngularMomentum() );
      CHECK( 0.0 == numbers.spin() );
      CHECK( 1.0 == numbers.totalAngularMomentum() );
      CHECK( +1 == numbers.parity() );

      auto radii = capture.radii();
      CHECK( 0.0 == Approx( radii.penetrabilityRadius( energy ).value ) );
      CHECK( 0.0 == Approx( radii.shiftFactorRadius( energy ).value ) );
      CHECK( 0.0 == Approx( radii.phaseShiftRadius( energy ).value ) );

      CHECK( 0.0 == capture.boundaryCondition() );

      CHECK( 0.0 == Approx( capture.Q().value ) );

      CHECK( 1.0 == Approx( capture.statisticalSpinFactor() ) );

      CHECK( false == Approx( capture.belowThreshold( energy ) ) );

      CHECK( 0.0 == Approx( capture.sommerfeldParameter( energy ) ) );
      CHECK( 0.0 == Approx( capture.waveNumber( energy ).value ) );
      CHECK( 1.0 == Approx( capture.penetrability( energy ) ) );
      CHECK( 0.0 == Approx( capture.shiftFactor( energy ) ) );
      CHECK( 0.0 == Approx( capture.phaseShift( energy ) ) );
      CHECK( 0.0 == Approx( capture.coulombPhaseShift( energy ) ) );
    } // THEN
  } // GIVEN
} // SCENARIO
