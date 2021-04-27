# standard imports
import unittest

# third party imports

# local imports
from elementary import ParticleID, ReactionID
from resonanceReconstruction import Particle, ParticlePair, ChannelQuantumNumbers, ChannelRadii
from resonanceReconstruction import PhotonChannel, NeutronChannel, ChargedParticleChannel, FissionChannel

class Test_r2_ParticleChannel( unittest.TestCase ) :
    """Unit test for the different channel classes."""

    def test_component( self ) :

        # particles
        photon = Particle( ParticleID( 'g' ), 0.0, 1., +1)
        neutron = Particle( ParticleID( 'n' ), 1.00866491582, 0.5, +1)
        proton = Particle( ParticleID( 'p' ), 1.00727647, 0.5, +1)
        cl36 = Particle( ParticleID( 'Cl36_e0' ), 35.968306822, 0., +1)
        cl35 = Particle( ParticleID( 'Cl35_e0' ), 34.968852694, 1.5, +1)
        cl35_e1 = Particle( ParticleID( 'Cl35_e1' ), 34.968852694, 1.5, +1)
        s36 = Particle( ParticleID( 'S36_e0' ), 35.967080699, 1.5, +1)

        # particle pairs
        elasticPair = ParticlePair( neutron, cl35 )
        inelasticPair = ParticlePair( neutron, cl35_e1 )
        capturePair = ParticlePair( photon, cl36 )
        protonEmissionPair = ParticlePair( proton, s36 )

        # Q values
        elasticQ = 0.0
        inelasticQ = -1.219440e+6
        captureQ = 0.0
        protonEmissionQ = 6.152200e+5

        # quantum numbers
        elasticNumbers = ChannelQuantumNumbers( 0, 1, 1, +1 )
        inelasticNumbers = ChannelQuantumNumbers( 0, 1, 1, +1 )
        captureNumbers = ChannelQuantumNumbers( 0, 0, 1, +1 )
        protonEmissionNumbers = ChannelQuantumNumbers( 0, 1, 1, +1 )

        # channel radii
        elasticRadii = ChannelRadii( 4.822220e-1, 3.667980e-1 )
        inelasticRadii = ChannelRadii( 4.822220e-1, 3.667980e-1 )
        captureRadii = ChannelRadii( 0.0 )
        protonEmissionRadii = ChannelRadii( 4.822220e-1, 3.667980e-1 )

        # custom channel identifier
        id = '1'

        # function inputs
        energy = 1e-5

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # photon channel, not an incident channel, no threshold
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        capture = PhotonChannel( incident = elasticPair, pair = capturePair,
                                 q_value = captureQ, numbers = captureNumbers,
                                 radii = captureRadii )

        self.assertEqual( 'photon,Cl36{0,0,1+}', capture.channel_id )
        self.assertEqual( 'n,Cl35->photon,Cl36', capture.reaction_id.symbol )

        self.assertEqual( False, capture.is_incident_channel )

        incident = capture.incident_particle_pair
        self.assertAlmostEqual( 1.00866491582, incident.particle.mass )
        self.assertAlmostEqual( 0.0, incident.particle.charge )
        self.assertEqual( 0.5, incident.particle.spin )
        self.assertEqual( +1, incident.particle.parity )
        self.assertAlmostEqual( 34.968852694, incident.residual.mass )
        self.assertAlmostEqual( 17. * 1.6021766208e-19, incident.residual.charge )
        self.assertEqual( 1.5, incident.residual.spin )
        self.assertEqual( +1, incident.residual.parity )
        self.assertEqual( 'n,Cl35', incident.pair_id.symbol )

        pair = capture.particle_pair
        self.assertAlmostEqual( 0.0, pair.particle.mass )
        self.assertAlmostEqual( 0.0, pair.particle.charge )
        self.assertEqual( 1.0, pair.particle.spin )
        self.assertEqual( +1, pair.particle.parity )
        self.assertAlmostEqual( 35.968306822, pair.residual.mass )
        self.assertAlmostEqual( 17. * 1.6021766208e-19, pair.residual.charge )
        self.assertEqual( 0.0, pair.residual.spin )
        self.assertEqual( +1, pair.residual.parity )
        self.assertEqual( 'photon,Cl36', pair.pair_id.symbol )

        numbers = capture.quantum_numbers
        self.assertEqual( 0, numbers.orbital_angular_momentum )
        self.assertEqual( 0.0, numbers.spin )
        self.assertEqual( 1.0, numbers.total_angular_momentum )
        self.assertEqual( +1, numbers.parity )

        radii = capture.radii
        self.assertAlmostEqual( 0.0, radii.penetrability_radius( energy ) )
        self.assertAlmostEqual( 0.0, radii.shift_factor_radius( energy ) )
        self.assertAlmostEqual( 0.0, radii.phase_shift_radius( energy ) )

        self.assertEqual( 0.0, capture.boundary_condition )

        self.assertAlmostEqual( 0.0, capture.Q )

        self.assertEqual( 1.0, capture.statistical_spin_factor )

        self.assertEqual( False, capture.below_threshold( energy ) )

        self.assertEqual( 0.0, capture.sommerfeld_parameter( energy ) )
        self.assertAlmostEqual( 0.0, capture.wave_number( energy ) )
        self.assertAlmostEqual( 1.0, capture.penetrability( energy ) )
        self.assertAlmostEqual( 0.0, capture.shift_factor( energy ) )
        self.assertAlmostEqual( 0.0, capture.phase_shift( energy ) )
        self.assertAlmostEqual( 0.0, capture.coulomb_phase_shift( energy ) )

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # neutron channel, incident channel, no threshold
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        elastic = NeutronChannel( incident = elasticPair, pair = elasticPair,
                                  q_value = elasticQ, numbers = elasticNumbers,
                                  radii = elasticRadii )

        self.assertEqual( 'n,Cl35{0,1,1+}', elastic.channel_id )
        self.assertEqual( 'n,Cl35->n,Cl35', elastic.reaction_id.symbol )

        self.assertEqual( True, elastic.is_incident_channel )

        incident = elastic.incident_particle_pair
        self.assertAlmostEqual( 1.00866491582, incident.particle.mass )
        self.assertAlmostEqual( 0.0, incident.particle.charge )
        self.assertEqual( 0.5, incident.particle.spin )
        self.assertEqual( +1, incident.particle.parity )
        self.assertAlmostEqual( 34.968852694, incident.residual.mass )
        self.assertAlmostEqual( 17. * 1.6021766208e-19, incident.residual.charge )
        self.assertEqual( 1.5, incident.residual.spin )
        self.assertEqual( +1, incident.residual.parity )
        self.assertEqual( 'n,Cl35', incident.pair_id.symbol )

        pair = elastic.particle_pair
        self.assertAlmostEqual( 1.00866491582, pair.particle.mass )
        self.assertAlmostEqual( 0.0, pair.particle.charge )
        self.assertEqual( 0.5, pair.particle.spin )
        self.assertEqual( +1, pair.particle.parity )
        self.assertAlmostEqual( 34.968852694, pair.residual.mass )
        self.assertAlmostEqual( 17. * 1.6021766208e-19, pair.residual.charge )
        self.assertEqual( 1.5, pair.residual.spin )
        self.assertEqual( +1, pair.residual.parity )
        self.assertEqual( 'n,Cl35', pair.pair_id.symbol )

        numbers = elastic.quantum_numbers
        self.assertEqual( 0, numbers.orbital_angular_momentum )
        self.assertEqual( 1.0, numbers.spin )
        self.assertEqual( 1.0, numbers.total_angular_momentum )
        self.assertEqual( +1, numbers.parity )

        radii = elastic.radii
        self.assertAlmostEqual( .4822220, radii.penetrability_radius( energy ) )
        self.assertAlmostEqual( .4822220, radii.shift_factor_radius( energy ) )
        self.assertAlmostEqual( .3667980, radii.phase_shift_radius( energy ) )

        self.assertEqual( 0.0, elastic.boundary_condition )

        self.assertAlmostEqual( 0.0, elastic.Q )

        self.assertEqual( 0.375, elastic.statistical_spin_factor )

        self.assertEqual( False, elastic.below_threshold( energy ) )

        self.assertAlmostEqual( 0.0, elastic.sommerfeld_parameter( energy ) )
        self.assertAlmostEqual( 6.75215238E-06, elastic.wave_number( energy ) )
        self.assertAlmostEqual( 3.25603642E-06, elastic.penetrability( energy ) )
        self.assertAlmostEqual( 0.0, elastic.shift_factor( energy ) )
        self.assertAlmostEqual( 2.47667599E-06, elastic.phase_shift( energy ) )
        self.assertAlmostEqual( 0.0, elastic.coulomb_phase_shift( energy ) )

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # neutron channel, not incident channel, threshold
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        inelastic = NeutronChannel( incident = elasticPair, pair = inelasticPair,
                                    q_value = inelasticQ, numbers = inelasticNumbers,
                                    radii = inelasticRadii )

        self.assertEqual( 'n,Cl35_e1{0,1,1+}', inelastic.channel_id )
        self.assertEqual( 'n,Cl35->n,Cl35_e1', inelastic.reaction_id.symbol )

        self.assertEqual( False, inelastic.is_incident_channel )

        incident = inelastic.incident_particle_pair
        self.assertAlmostEqual( 1.00866491582, incident.particle.mass )
        self.assertAlmostEqual( 0.0, incident.particle.charge )
        self.assertEqual( 0.5, incident.particle.spin )
        self.assertEqual( +1, incident.particle.parity )
        self.assertAlmostEqual( 34.968852694, incident.residual.mass )
        self.assertAlmostEqual( 17. * 1.6021766208e-19, incident.residual.charge )
        self.assertEqual( 1.5, incident.residual.spin )
        self.assertEqual( +1, incident.residual.parity )
        self.assertEqual( 'n,Cl35', incident.pair_id.symbol )

        pair = inelastic.particle_pair
        self.assertAlmostEqual( 1.00866491582, pair.particle.mass )
        self.assertAlmostEqual( 0.0, pair.particle.charge )
        self.assertEqual( 0.5, pair.particle.spin )
        self.assertEqual( +1, pair.particle.parity )
        self.assertAlmostEqual( 34.968852694, pair.residual.mass )
        self.assertAlmostEqual( 17. * 1.6021766208e-19, pair.residual.charge )
        self.assertEqual( 1.5, pair.residual.spin )
        self.assertEqual( +1, pair.residual.parity )
        self.assertEqual( 'n,Cl35_e1', pair.pair_id.symbol )

        numbers = inelastic.quantum_numbers
        self.assertEqual( 0, numbers.orbital_angular_momentum )
        self.assertEqual( 1.0, numbers.spin )
        self.assertEqual( 1.0, numbers.total_angular_momentum )
        self.assertEqual( +1, numbers.parity )

        radii = inelastic.radii
        self.assertAlmostEqual( .4822220, radii.penetrability_radius( energy ) )
        self.assertAlmostEqual( .4822220, radii.shift_factor_radius( energy ) )
        self.assertAlmostEqual( .3667980, radii.phase_shift_radius( energy ) )

        self.assertEqual( 0.0, inelastic.boundary_condition )

        self.assertAlmostEqual( -1.219440e+6, inelastic.Q )

        self.assertEqual( 0.375, inelastic.statistical_spin_factor )

        self.assertEqual( True, inelastic.below_threshold( energy ) )
        self.assertEqual( False, inelastic.below_threshold( 15e+6 ) )

        self.assertAlmostEqual( 0.0, inelastic.sommerfeld_parameter( energy ) )
        self.assertAlmostEqual( 2.39164854E+00, inelastic.wave_number( energy ) )
        self.assertAlmostEqual( 1.15330554E+00, inelastic.penetrability( energy ) )
        self.assertAlmostEqual( 0.0, inelastic.shift_factor( energy ) )
        self.assertAlmostEqual( 8.77251899E-01, inelastic.phase_shift( energy ) )
        self.assertAlmostEqual( 0.0, inelastic.coulomb_phase_shift( energy ) )

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # proton channel, not an incident channel, no threshold
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        protonEmission = ChargedParticleChannel(
                             incident = elasticPair, pair = protonEmissionPair,
                             q_value = protonEmissionQ, numbers = protonEmissionNumbers,
                             radii = protonEmissionRadii )

        self.assertEqual( 'p,S36{0,1,1+}', protonEmission.channel_id )
        self.assertEqual( 'n,Cl35->p,S36', protonEmission.reaction_id.symbol )

        self.assertEqual( False, protonEmission.is_incident_channel )

        incident = protonEmission.incident_particle_pair
        self.assertAlmostEqual( 1.00866491582, incident.particle.mass )
        self.assertAlmostEqual( 0.0, incident.particle.charge )
        self.assertEqual( 0.5, incident.particle.spin )
        self.assertEqual( +1, incident.particle.parity )
        self.assertAlmostEqual( 34.968852694, incident.residual.mass )
        self.assertAlmostEqual( 17. * 1.6021766208e-19, incident.residual.charge )
        self.assertEqual( 1.5, incident.residual.spin )
        self.assertEqual( +1, incident.residual.parity )
        self.assertEqual( 'n,Cl35', incident.pair_id.symbol )

        pair = protonEmission.particle_pair
        self.assertAlmostEqual( 1.00727647, pair.particle.mass )
        self.assertAlmostEqual( 1.6021766208e-19, pair.particle.charge )
        self.assertEqual( 0.5, pair.particle.spin )
        self.assertEqual( +1, pair.particle.parity )
        self.assertAlmostEqual( 35.967080699, pair.residual.mass )
        self.assertAlmostEqual( 16. * 1.6021766208e-19, pair.residual.charge )
        self.assertEqual( 1.5, pair.residual.spin )
        self.assertEqual( +1, pair.residual.parity )
        self.assertEqual( 'p,S36', pair.pair_id.symbol )

        numbers = protonEmission.quantum_numbers
        self.assertEqual( 0, numbers.orbital_angular_momentum )
        self.assertEqual( 1.0, numbers.spin )
        self.assertEqual( 1.0, numbers.total_angular_momentum )
        self.assertEqual( +1, numbers.parity )

        radii = protonEmission.radii
        self.assertAlmostEqual( .4822220, radii.penetrability_radius( energy ) )
        self.assertAlmostEqual( .4822220, radii.shift_factor_radius( energy ) )
        self.assertAlmostEqual( .3667980, radii.phase_shift_radius( energy ) )

        self.assertEqual( 0.0, protonEmission.boundary_condition )

        self.assertAlmostEqual( 6.152200e+5, protonEmission.Q )

        self.assertEqual( 0.375, protonEmission.statistical_spin_factor )

        self.assertEqual( False, protonEmission.below_threshold( energy ) )

        self.assertAlmostEqual( 3.17996084E+00, protonEmission.sommerfeld_parameter( energy ) )
        self.assertAlmostEqual( 1.69828445E+00, protonEmission.wave_number( energy ) )
        self.assertAlmostEqual( 0.000027793, protonEmission.penetrability( energy ) )
        self.assertAlmostEqual( -1.8734549658, protonEmission.shift_factor( energy ) )
        self.assertAlmostEqual( 0.0000020468, protonEmission.phase_shift( energy ) )
        self.assertAlmostEqual( 0.0, protonEmission.coulomb_phase_shift( energy ) )

if __name__ == '__main__' :

    unittest.main()
