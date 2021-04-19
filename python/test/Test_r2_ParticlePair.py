# standard imports
import unittest

# third party imports

# local imports
from elementary import ParticleID, ParticlePairID
from resonanceReconstruction import Particle, ParticlePair

class Test_r2_ParticlePair( unittest.TestCase ) :
    """Unit test for the ParticlePair class."""

    def test_component( self ) :

        # particles
        photon = Particle( ParticleID( 'g' ), 0.0, 1., +1)
        neutron = Particle( ParticleID( 'n' ), 1.00866491582, 0.5, +1)
        proton = Particle( ParticleID( 'p' ), 1.00727647, 0.5, +1)
        cl36 = Particle( ParticleID( 'Cl36_e0' ), 35.968306822, 0., +1)
        cl35 = Particle( ParticleID( 'Cl35_e0' ), 34.968852694, 1.5, +1)
        s36 = Particle( ParticleID( 'S36_e0' ), 35.967080699, 1.5, +1)

        # custom particle pair identifier
        id = ParticlePairID( 'fission' );

        # a ParticlePair can be constructed without a pair ID
        pair = ParticlePair( particle = neutron, residual = cl35 )

        self.assertAlmostEqual( 1.00866491582, pair.particle.mass )
        self.assertAlmostEqual( 0.0, pair.particle.charge )
        self.assertEqual( 0.5, pair.particle.spin )
        self.assertEqual( +1, pair.particle.parity )

        self.assertAlmostEqual( 34.968852694, pair.residual.mass )
        self.assertAlmostEqual( 17.0 * 1.60217662e-19, pair.residual.charge )
        self.assertEqual( 1.5, pair.residual.spin )
        self.assertEqual( +1, pair.residual.parity )

        self.assertEqual( 'n,Cl35', pair.pair_id.symbol )

        self.assertAlmostEqual( 34.968852694 / ( 34.968852694 + 1.00866491582 ),
                                pair.mass_ratio )
        self.assertEqual( 1.00866491582 * 34.968852694 / ( 34.968852694 + 1.00866491582 ),
                          pair.reduced_mass )

        # a ParticlePair can be constructed with a pair ID
        pair = ParticlePair( neutron, cl35, id );

        self.assertAlmostEqual( 1.00866491582, pair.particle.mass )
        self.assertAlmostEqual( 0.0, pair.particle.charge )
        self.assertEqual( 0.5, pair.particle.spin )
        self.assertEqual( +1, pair.particle.parity )

        self.assertAlmostEqual( 34.968852694, pair.residual.mass )
        self.assertAlmostEqual( 17.0 * 1.60217662e-19, pair.residual.charge )
        self.assertEqual( 1.5, pair.residual.spin )
        self.assertEqual( +1, pair.residual.parity )

        self.assertEqual( 'fission', pair.pair_id.symbol )

        self.assertEqual( 34.968852694 / ( 34.968852694 + 1.00866491582 ),
                          pair.mass_ratio )
        self.assertEqual( 1.00866491582 * 34.968852694 / ( 34.968852694 + 1.00866491582 ),
                          pair.reduced_mass )

if __name__ == '__main__' :

    unittest.main()
