import os
from ase.clease import CEBulk, NewStructures, Concentration
from ase.build import bulk
from ase.clease.tools import wrap_and_sort_by_position

db_name = 'fu.db'

def test_1():
    basis_elements = [['Li', 'Ru', 'X'], ['O', 'X']]
    A_eq = [[0, 3, 0, 0, 0]]
    b_eq = [1]
    A_lb = [[0, 0, 0, 3, 0]]
    b_lb = [2]
    conc = Concentration(basis_elements=basis_elements, A_eq=A_eq, b_eq=b_eq,
                         A_lb=A_lb, b_lb=b_lb)
    setting = CEBulk(crystalstructure="rocksalt",
                     a=4.0,
                     size=[2, 2, 2],
                     concentration=conc,
                     db_name=db_name,
                     max_cluster_size=3,
                     max_cluster_dia=4.)
    newstruct = NewStructures(setting=setting)

    test = bulk(name='LiO', crystalstructure='rocksalt', a=4.0) * (2,2,2)
    atoms = wrap_and_sort_by_position(test.copy())
    atoms[0].symbol = 'Ru'
    atoms[1].symbol = 'X'
    atoms[14].symbol = 'X'
    fu = newstruct._get_formula_unit(atoms)
    assert fu == "Li6Ru1X1_O7X1"

    atoms = wrap_and_sort_by_position(test.copy())
    fu = newstruct._get_formula_unit(atoms)
    assert fu == "Li1_O1"

    atoms = wrap_and_sort_by_position(test.copy())
    replace_cat = [0, 1, 2, 5]
    replace_an = [4, 7, 9, 12]
    for cat in replace_cat:
        atoms[cat].symbol = 'X'
    for an in replace_an:
        atoms[an].symbol = 'X'

    fu = newstruct._get_formula_unit(atoms)
    assert fu == "Li1X1_O1X1"
    os.remove(db_name)

 

test_1()
