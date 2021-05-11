# standard imports
import unittest

# third party imports

# local imports
from resonanceReconstruction import Boundary

class Test_r2_Boundary( unittest.TestCase ) :
    """Unit test for the Boundary enumerator."""

    def test_component( self ) :

        option = Boundary.ShiftFactor

        self.assertTrue( option == Boundary.ShiftFactor )
        self.assertFalse( option == Boundary.Constant )

        self.assertFalse( option != Boundary.ShiftFactor )
        self.assertTrue( option != Boundary.Constant )

        option = Boundary.Constant

        self.assertFalse( option == Boundary.ShiftFactor )
        self.assertTrue( option == Boundary.Constant )

        self.assertTrue( option != Boundary.ShiftFactor )
        self.assertFalse( option != Boundary.Constant )

if __name__ == '__main__' :

    unittest.main()
