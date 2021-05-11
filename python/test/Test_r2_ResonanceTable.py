# standard imports
import unittest

# third party imports

# local imports
from resonanceReconstruction import Resonance, ResonanceTable

class Test_r2_Resonance( unittest.TestCase ) :
    """Unit test for the Resonance class."""

    def test_component( self ) :

        # resonance data
        channels = [ "1", "2" ]
        resonances = [ Resonance( 6.823616e+4, [ 2.179040e+2, 1.000000e-5 ], 3.933600e-1 ),
                       Resonance( 1.150980e+5, [ 4.307780e+0, 0.0 ], 7.390000e-1 ),
                       Resonance( 1.825230e+5, [ 1.759740e+3, 4.000000e-1 ], 7.451500e-1 ) ]

        table = ResonanceTable( channels = channels, resonances = resonances )

        self.assertEqual( 2, table.number_channels )
        self.assertEqual( 2, len( table.channels ) )
        self.assertEqual( "1", table.channels[0] );
        self.assertEqual( "2", table.channels[1] );

        self.assertEqual( 3, table.number_resonances )
        self.assertEqual( 3, len( table.resonances ) )
        self.assertEqual( 3, len( table.energies ) )
        self.assertAlmostEqual( 6.823616e+4, table.energies[0] )
        self.assertAlmostEqual( 1.150980e+5, table.energies[1] )
        self.assertAlmostEqual( 1.825230e+5, table.energies[2] )

        resonance = table.resonances[0]
        self.assertAlmostEqual( 6.823616e+4, resonance.energy )
        self.assertAlmostEqual( 3.933600e-1, resonance.eliminated_width )
        self.assertEqual( 2, len( resonance.widths ) )
        self.assertAlmostEqual( 2.179040e+2, resonance.widths[0] )
        self.assertAlmostEqual( 1.000000e-5, resonance.widths[1] )

        resonance = table.resonances[1]
        self.assertAlmostEqual( 1.150980e+5, resonance.energy )
        self.assertAlmostEqual( 7.390000e-1, resonance.eliminated_width )
        self.assertEqual( 2, len( resonance.widths ) )
        self.assertAlmostEqual( 4.307780e+0, resonance.widths[0] )
        self.assertAlmostEqual( 0.0, resonance.widths[1] )

        resonance = table.resonances[2]
        self.assertAlmostEqual( 1.825230e+5, resonance.energy )
        self.assertAlmostEqual( 7.451500e-1, resonance.eliminated_width )
        self.assertEqual( 2, len( resonance.widths ) )
        self.assertAlmostEqual( 1.759740e+3, resonance.widths[0] )
        self.assertAlmostEqual( 4.000000e-1, resonance.widths[1] )

    def test_failures( self ) :

        print( '\n' )

        # wrong number of widths in a Resonance
        with self.assertRaises( Exception ) :

            channels = [ "1", "2" ]
            resonances = [ Resonance( 6.823616e+4, [ 2.179040e+2, 1.000000e-5 ], 3.933600e-1 ),
                           Resonance( 1.150980e+5, [ 4.307780e+0 ], 7.390000e-1 ),
                           Resonance( 1.825230e+5, [ 1.759740e+3, 4.000000e-1 ], 7.451500e-1 ) ]

            table = ResonanceTable( channels = channels, resonances = resonances )

if __name__ == '__main__' :

    unittest.main()
