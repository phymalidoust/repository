import unittest
#from ase import Atoms
from ase.build import molecule
from ase.calculators import lj

class MonteCarloTest( unittest.TestCase ):
    def test_no_except( self ):
        try:
            from ase.montecarlo.montecarlo import Montecarlo
            # Create an atoms objec consisting of more than 2 atoms
            atoms = molecule( "H2O" )

            # Use the LennardJones calculator for simplicity
            calc = lj.LennardJones()
            atoms.set_calculator( calc )

            mc = Montecarlo( atoms, 200.0 )

            # Run 10 steps to verify the code runs without throwing exceptions
            mc.runMC( steps=10, verbose=False )
        except Exception as exc:
            self.fail( str(exc) )

if __name__ == "__main__":
    unittest.main()
