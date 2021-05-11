# standard imports
import unittest

# third party imports

# local imports
from resonanceReconstruction import Resonance

class Test_r2_Resonance( unittest.TestCase ) :
    """Unit test for the Resonance class."""

    def test_component( self ) :

        # single resonance data
        energy = 6.823616e+4
        eliminated = 3.933600e-1
        widths = [ 2.179040e+2, 1.000000e-5 ]

        # a Resonance can be constructed with an eliminated width

        resonance = Resonance( energy = energy, widths = widths,
                               eliminated = eliminated )

        self.assertAlmostEqual( 6.823616e+4, resonance.energy )

        self.assertAlmostEqual( 3.933600e-1, resonance.eliminated_width )

        self.assertEqual( 2, len( resonance.widths ) )
        self.assertAlmostEqual( 2.179040e+2, resonance.widths[0] )
        self.assertAlmostEqual( 1.000000e-5, resonance.widths[1] )

        # a Resonance can be constructed without an eliminated width

        resonance = Resonance( energy = energy, widths = widths )

        self.assertAlmostEqual( 6.823616e+4, resonance.energy )

        self.assertEqual( 0.0, resonance.eliminated_width )

        self.assertEqual( 2, len( resonance.widths ) )
        self.assertAlmostEqual( 2.179040e+2, resonance.widths[0] )
        self.assertAlmostEqual( 1.000000e-5, resonance.widths[1] )

if __name__ == '__main__' :

    unittest.main()
