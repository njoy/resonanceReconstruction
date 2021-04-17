# standard imports
import unittest

# third party imports

# local imports
from elementary import ParticleID
from resonanceReconstruction import Particle

class Test_r2_Particle( unittest.TestCase ) :
    """Unit test for the Particle class."""

    def test_component( self ) :

        # neutron
        neutron = Particle( id = ParticleID( "n" ), mass = 1.008664,
                            charge = 0.0, spin = 0.5, parity = +1 )

        # proton
        proton = Particle( id = ParticleID( "p" ), mass = 1.007276,
                           charge = 1.60217662e-19, spin = 0.5, parity = +1 )

        # U235
        u235 = Particle( id = ParticleID( "U235" ), mass = 235.0439299,
                         charge = 0.0, spin = 0., parity = +1 )

        self.assertEqual( ParticleID( "n" ), neutron.particle_id )
        self.assertEqual( 1.008664, neutron.mass )
        self.assertEqual( 0., neutron.charge )
        self.assertEqual( 0.5, neutron.spin )
        self.assertEqual( +1, neutron.parity )

        self.assertEqual( ParticleID( "p" ), proton.particle_id )
        self.assertEqual( 1.007276, proton.mass )
        self.assertEqual( 1.60217662e-19, proton.charge )
        self.assertEqual( 0.5, proton.spin )
        self.assertEqual( +1, proton.parity )

        self.assertEqual( ParticleID( "U235" ), u235.particle_id )
        self.assertEqual( 235.0439299, u235.mass )
        self.assertEqual( 0., u235.charge )
        self.assertEqual( 0., u235.spin )
        self.assertEqual( +1, u235.parity )

    def test_failures( self ) :

        print( '\n' )

        # negative mass value
        with self.assertRaises( Exception ) :

            neutron = Particle( id = ParticleID( "n" ), mass = -1.008664,
                                charge = 0.0, spin = 0.5, parity = +1 )

        # negative charge value
        with self.assertRaises( Exception ) :

            neutron = Particle( id = ParticleID( "n" ), mass = 1.008664,
                                charge = -1.0, spin = 0.5, parity = +1 )

if __name__ == '__main__' :

    unittest.main()
