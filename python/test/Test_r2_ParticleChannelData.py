# standard imports
import unittest

# third party imports

# local imports
from elementary import ParticleID, ReactionID
from resonanceReconstruction import Particle, ParticlePair, ChannelQuantumNumbers, ChannelRadii
from resonanceReconstruction import PhotonChannel, NeutronChannel, ChargedParticleChannel, FissionChannel
from resonanceReconstruction import ParticleChannelData

class Test_r2_ParticleChannelData( unittest.TestCase ) :
    """Unit test for the ParticleChannelData classes."""

    def test_component( self ) :

        # particles
        neutron = Particle( ParticleID( 'n' ), 1.00866491582, 0.5, +1 )
        cl35 = Particle( ParticleID( 'Cl35_e0' ), 34.968852694, 1.5, +1 )

        # particle pair
        pair = ParticlePair( neutron, cl35 )

        # quantum numbers
        numbers = ChannelQuantumNumbers( 0, 1.0, 1.0, +1 )

        # channel radii
        radii = ChannelRadii( 4.822220e-1, 3.667980e-1 )

        # channel
        elastic = NeutronChannel( incident = pair, pair = pair,
                                  q_value = 0., numbers = numbers,
                                  radii = radii )

        # width and energy data
        energies = [ 1.2345e+1, 3500. ]
        widths = [ 3.4e-2, 6.5e-1 ]

        data = ParticleChannelData( channel = elastic,
                                    energies = energies, widths = widths )

        self.assertEqual( 'n,Cl35{0,1,1+}', data.channel_id )
        self.assertEqual( 'n,Cl35->n,Cl35', data.reaction_id.symbol )

        self.assertEqual( True, data.is_incident_channel )
        self.assertEqual( False, data.is_eliminated_channel )

        incident = data.incident_particle_pair
        self.assertAlmostEqual( 1.00866491582, incident.particle.mass )
        self.assertAlmostEqual( 0.0, incident.particle.charge )
        self.assertAlmostEqual( 0.5, incident.particle.spin )
        self.assertEqual( +1, incident.particle.parity )
        self.assertAlmostEqual( 34.968852694, incident.residual.mass )
        self.assertAlmostEqual( 17. * 1.6021766208e-19, incident.residual.charge )
        self.assertAlmostEqual( 1.5, incident.residual.spin )
        self.assertEqual( +1, incident.residual.parity )
        self.assertEqual( 'n,Cl35', incident.pair_id.symbol )

        pair = data.particle_pair
        self.assertAlmostEqual( 1.00866491582, pair.particle.mass )
        self.assertAlmostEqual( 0.0, pair.particle.charge )
        self.assertAlmostEqual( 0.5, pair.particle.spin )
        self.assertEqual( +1, pair.particle.parity )
        self.assertAlmostEqual( 34.968852694, pair.residual.mass )
        self.assertAlmostEqual( 17. * 1.6021766208e-19, pair.residual.charge )
        self.assertAlmostEqual( 1.5, pair.residual.spin )
        self.assertEqual( +1, pair.residual.parity )
        self.assertEqual( 'n,Cl35', pair.pair_id.symbol )

        numbers = data.quantum_numbers
        self.assertEqual( 0, numbers.orbital_angular_momentum )
        self.assertEqual( 1.0, numbers.spin )
        self.assertEqual( 1.0, numbers.total_angular_momentum )
        self.assertEqual( +1, numbers.parity )

        radii = data.radii
        energy = 1e-5
        self.assertAlmostEqual( .4822220, radii.penetrability_radius( energy ) )
        self.assertAlmostEqual( .4822220, radii.shift_factor_radius( energy ) )
        self.assertAlmostEqual( .3667980, radii.phase_shift_radius( energy ) )

        self.assertEqual( 0.0, data.boundary_condition )

        self.assertEqual( 0.0, data.Q )

        self.assertEqual( 2, len( data.energies ) )
        self.assertAlmostEqual( 12.345, data.energies[0] )
        self.assertAlmostEqual( 3500., data.energies[1] )
        self.assertEqual( 2, len( data.widths ) )
        self.assertAlmostEqual( 3.4e-2, data.widths[0] )
        self.assertAlmostEqual( 0.65, data.widths[1] )

if __name__ == '__main__' :

    unittest.main()
