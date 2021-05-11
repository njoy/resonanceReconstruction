# standard imports
import unittest

# third party imports

# local imports
from resonanceReconstruction import Formalism

class Test_r2_Formalism( unittest.TestCase ) :
    """Unit test for the Formalism enumerator."""

    def test_component( self ) :

        option = Formalism.SingleLevelBreitWigner

        self.assertTrue( option == Formalism.SingleLevelBreitWigner )
        self.assertFalse( option == Formalism.MultiLevelBreitWigner )
        self.assertFalse( option == Formalism.ReichMoore )
        self.assertFalse( option == Formalism.GeneralRMatrix )

        self.assertFalse( option != Formalism.SingleLevelBreitWigner )
        self.assertTrue( option != Formalism.MultiLevelBreitWigner )
        self.assertTrue( option != Formalism.ReichMoore )
        self.assertTrue( option != Formalism.GeneralRMatrix )

        option = Formalism.MultiLevelBreitWigner

        self.assertFalse( option == Formalism.SingleLevelBreitWigner )
        self.assertTrue( option == Formalism.MultiLevelBreitWigner )
        self.assertFalse( option == Formalism.ReichMoore )
        self.assertFalse( option == Formalism.GeneralRMatrix )

        self.assertTrue( option != Formalism.SingleLevelBreitWigner )
        self.assertFalse( option != Formalism.MultiLevelBreitWigner )
        self.assertTrue( option != Formalism.ReichMoore )
        self.assertTrue( option != Formalism.GeneralRMatrix )

        option = Formalism.ReichMoore

        self.assertFalse( option == Formalism.SingleLevelBreitWigner )
        self.assertFalse( option == Formalism.MultiLevelBreitWigner )
        self.assertTrue( option == Formalism.ReichMoore )
        self.assertFalse( option == Formalism.GeneralRMatrix )

        self.assertTrue( option != Formalism.SingleLevelBreitWigner )
        self.assertTrue( option != Formalism.MultiLevelBreitWigner )
        self.assertFalse( option != Formalism.ReichMoore )
        self.assertTrue( option != Formalism.GeneralRMatrix )

        option = Formalism.GeneralRMatrix

        self.assertFalse( option == Formalism.SingleLevelBreitWigner )
        self.assertFalse( option == Formalism.MultiLevelBreitWigner )
        self.assertFalse( option == Formalism.ReichMoore )
        self.assertTrue( option == Formalism.GeneralRMatrix )

        self.assertTrue( option != Formalism.SingleLevelBreitWigner )
        self.assertTrue( option != Formalism.MultiLevelBreitWigner )
        self.assertTrue( option != Formalism.ReichMoore )
        self.assertFalse( option != Formalism.GeneralRMatrix )

if __name__ == '__main__' :

    unittest.main()
