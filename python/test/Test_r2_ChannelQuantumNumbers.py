# standard imports
import unittest

# third party imports

# local imports
from resonanceReconstruction import ChannelQuantumNumbers

class Test_r2_ChannelQuantumNumbers( unittest.TestCase ) :
    """Unit test for the ChannelQuantumNumbers class."""

    def test_component( self ) :

        numbers = ChannelQuantumNumbers( l = 1, s = 0.5, J = 1.5, parity = +1 )

        self.assertEqual( 1, numbers.orbital_angular_momentum )
        self.assertEqual( 0.5, numbers.spin )
        self.assertEqual( 1.5, numbers.total_angular_momentum )
        self.assertEqual( +1, numbers.parity )
        self.assertEqual( "{1,1/2,3/2+}", str( numbers ) )

if __name__ == '__main__' :

    unittest.main()
