# standard imports
import unittest

# third party imports

# local imports
from resonanceReconstruction import ChannelRadii

class Test_r2_ChannelRadii( unittest.TestCase ) :
    """Unit test for the ChannelRadiir class."""

    def test_component( self ) :

        # radii
        radius = 0.1

        trueRadius = 0.2
        effectiveRadius = 0.3

        pRadius = 0.4
        sRadius = 0.5
        phiRadius = 0.6

        # energy
        energy = 1e-5;

        radii = ChannelRadii( radius = radius );

        self.assertAlmostEqual( 0.1, radii.penetrability_radius( energy ) )
        self.assertAlmostEqual( 0.1, radii.shift_factor_radius( energy ) )
        self.assertAlmostEqual( 0.1, radii.phase_shift_radius( energy ) )

        radii = ChannelRadii( true_radius = trueRadius,
                              effective_radius = effectiveRadius );

        self.assertAlmostEqual( 0.2, radii.penetrability_radius( energy ) )
        self.assertAlmostEqual( 0.2, radii.shift_factor_radius( energy ) )
        self.assertAlmostEqual( 0.3, radii.phase_shift_radius( energy ) )

        radii = ChannelRadii( penetrability_radius = pRadius,
                              shift_factor_radius = sRadius,
                              phase_shift_radius = phiRadius );

        self.assertAlmostEqual( 0.4, radii.penetrability_radius( energy ) )
        self.assertAlmostEqual( 0.5, radii.shift_factor_radius( energy ) )
        self.assertAlmostEqual( 0.6, radii.phase_shift_radius( energy ) )

if __name__ == '__main__' :

    unittest.main()
