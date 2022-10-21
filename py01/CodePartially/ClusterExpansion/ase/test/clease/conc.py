import numpy as np
from ase.clease.concentration import Concentration
from collections import OrderedDict


def test_full_range():
    basis_elements = [['Li', 'Ru', 'X'], ['O', 'X']]
    conc = Concentration(basis_elements=basis_elements)
    conc = conc.get_random_concentration()
    sum1 = np.sum(conc[:3])
    assert abs(sum1 - 1) < 1E-9
    sum2 = np.sum(conc[3:])
    assert abs(sum2 - 1) < 1E-9


def fixed_composition():
    basis_elements = [['Li', 'Ru'], ['O', 'X']]
    A_eq = [[0, 3, 0, 0], [0, 0, 0, 2]]
    b_eq = [1, 1]
    conc = Concentration(basis_elements=basis_elements, A_eq=A_eq, b_eq=b_eq)
    rand = conc.get_random_concentration()
    assert np.allclose(rand, np.array([2./3, 1./3, 0.5, 0.5]))


def test_fix_Ru_composition():
    basis_elements = [['Li', 'Ru', 'X'], ['O', 'X']]
    A_eq = [[0, 3, 0, 0, 0]]
    b_eq = [1]
    A_lb = [[0, 0, 0, 3, 0]]
    b_lb = [2]
    conc = Concentration(basis_elements=basis_elements, A_eq=A_eq, b_eq=b_eq,
                         A_lb=A_lb, b_lb=b_lb)
    rand = conc.get_random_concentration()
    sum1 = np.sum(rand[:3])
    assert abs(sum1 - 1) < 1E-9
    sum2 = np.sum(rand[3:])
    assert abs(sum2 - 1) < 1E-9


def test_conc_range():
    basis_elements = [['Li', 'Ru', 'X'], ['O', 'X']]
    conc = Concentration(basis_elements=basis_elements)
    ranges = [[(0, 1), (1./3, 1./3), (0, 1)], [(2./3, 1), (0, 1)]]
    conc.set_conc_ranges(ranges)
    rand = conc.get_random_concentration()
    sum1 = np.sum(rand[:3])
    assert abs(sum1 - 1) < 1E-9
    sum2 = np.sum(rand[3:])
    assert abs(sum2 - 1) < 1E-9
    assert abs(rand[1] - 1./3) < 1E-9
    assert rand[3] > 2./3


def test_formula_unit1():
    basis_elements = [['Li', 'Ru', 'X'], ['O', 'X']]
    conc = Concentration(basis_elements=basis_elements)
    # test _get_integers() function
    example_str = "df<12> kd542,32"
    expected_ints = [12, 542, 32]
    assert conc._get_integers(example_str) == expected_ints


def test_formula_unit2():
    basis_elements = [['Li', 'Ru', 'X'], ['O', 'X']]
    formulas = ["Li<x>Ru<1>X<2-x>", "O<3-y>X<y>"]
    var_range = OrderedDict()
    var_range["x"] = (0, 2)
    var_range["y"] = (0, 0.75)
    conc = Concentration(basis_elements=basis_elements)
    conc.set_conc_formula_unit(formulas=formulas, variable_range=var_range)
    A_eq = np.array([[1, 1, 1, 0, 0],
                     [0, 0, 0, 1, 1],
                     [0, 3, 0, 0, 0]])
    b_eq = np.array([1, 1, 1])
    assert np.allclose(A_eq, conc.A_eq)
    assert np.allclose(b_eq, conc.b_eq)
    A_lb_implicit = np.identity(5)
    b_lb_implicit = np.zeros(5)
    A_lb = np.array([[3, 0, 0, 0, 0],
                     [-3, 0, 0, 0, 0],
                     [0, 0, 0, 0, 3],
                     [0, 0, 0, 0, -3]])
    b_lb = np.array([0, -2, 0, -0.75])
    A_lb = np.vstack((A_lb_implicit, A_lb))
    b_lb = np.append(b_lb_implicit, b_lb)
    assert np.allclose(A_lb, conc.A_lb)
    assert np.allclose(b_lb, conc.b_lb)


def test_formula_unit3():
    basis_elements = [['Li', 'V', 'X'], ['O', 'F']]
    formulas = ["Li<x>V<1>X<2-x>", "O<2>F<1>"]
    range = OrderedDict({"x": (0, 2)})
    conc = Concentration(basis_elements=basis_elements)
    conc.set_conc_formula_unit(formulas=formulas, variable_range=range)
    A_eq = np.array([[1, 1, 1, 0, 0],
                     [0, 0, 0, 1, 1],
                     [0, 3, 0, 0, 0],
                     [0, 0, 0, 3, 0],
                     [0, 0, 0, 0, 3]])
    b_eq = np.array([1, 1, 1, 2, 1])
    assert np.allclose(A_eq, conc.A_eq)
    assert np.allclose(b_eq, conc.b_eq)
    A_lb_implicit = np.identity(5)
    b_lb_implicit = np.zeros(5)
    A_lb = np.array([[3, 0, 0, 0, 0],
                     [-3, 0, 0, 0, 0]])
    b_lb = np.array([0, -2])
    A_lb = np.vstack((A_lb_implicit, A_lb))
    b_lb = np.append(b_lb_implicit, b_lb)
    assert np.allclose(A_lb, conc.A_lb)
    assert np.allclose(b_lb, conc.b_lb)


def test_formula_unit4():
    basis_elements = [['Al', 'Mg', 'Si']]
    formulas = ['Al<4-4x>Mg<3x>Si<x>']
    range = OrderedDict({"x": (0, 1)})
    conc = Concentration(basis_elements=basis_elements)
    conc.set_conc_formula_unit(formulas=formulas, variable_range=range)
    A_eq = np.array([[1, 1, 1],
                     [0, 1, -3]])
    b_eq = np.array([1, 0])
    assert np.allclose(A_eq, conc.A_eq)
    assert np.allclose(b_eq, conc.b_eq)
    A_lb_implicit = np.identity(3)
    b_lb_implicit = np.zeros(3)
    A_lb = np.array([[0, 0, 4],
                     [0, 0, -4]])
    b_lb = np.array([0, -1])
    A_lb = np.vstack((A_lb_implicit, A_lb))
    b_lb = np.append(b_lb_implicit, b_lb)
    assert np.allclose(A_lb, conc.A_lb)
    assert np.allclose(b_lb, conc.b_lb)


def test_formula_unit5():
    basis_elements = [['Al', 'Mg', 'Si']]
    formulas = ['Al<2-2x>Mg<x>Si<x>']
    range = OrderedDict({"x": (0, 1)})
    conc = Concentration(basis_elements=basis_elements)
    conc.set_conc_formula_unit(formulas=formulas, variable_range=range)
    A_eq = np.array([[1, 1, 1],
                     [0, -1, 1]])
    b_eq = np.array([1, 0])
    assert np.allclose(A_eq, conc.A_eq)
    assert np.allclose(b_eq, conc.b_eq)
    A_lb_implicit = np.identity(3)
    b_lb_implicit = np.zeros(3)
    A_lb = np.array([[0, 2, 0],
                     [0, -2, 0]])
    b_lb = np.array([0, -1])
    A_lb = np.vstack((A_lb_implicit, A_lb))
    b_lb = np.append(b_lb_implicit, b_lb)
    assert np.allclose(A_lb, conc.A_lb)
    assert np.allclose(b_lb, conc.b_lb)


def test_formula_unit6():
    basis_elements = [['Al', 'Mg', 'Si', 'X']]
    formulas = ['Al<3-x-2y>Mg<y>Si<y>X<x>']
    range = {"x": (0, 1), "y": (0, 1)}
    conc = Concentration(basis_elements=basis_elements)
    conc.set_conc_formula_unit(formulas=formulas, variable_range=range)
    A_eq = np.array([[1, 1, 1, 1],
                     [0, -1, 1, 0]])
    b_eq = np.array([1, 0])
    assert np.allclose(A_eq, conc.A_eq)
    assert np.allclose(b_eq, conc.b_eq)
    A_lb_implicit = np.identity(4)
    b_lb_implicit = np.zeros(4)
    A_lb = np.array([[0, 0, 0, 3],
                     [0, 0, 0, -3],
                     [0, 3, 0, 0],
                     [0, -3, 0, 0]])
    b_lb = np.array([0, -1, 0, -1])
    A_lb = np.vstack((A_lb_implicit, A_lb))
    b_lb = np.append(b_lb_implicit, b_lb)
    assert np.allclose(A_lb, conc.A_lb)
    assert np.allclose(b_lb, conc.b_lb)

def test_three_interlinked_basis():
    basis_elements = [["Al", "Mg", "Si"], ["X", "O"], ["Ta", "Se"]]
    A_eq = [[1, 0, 0, -1, 0, 0, 0]]
    b_eq = [0]
    A_lb = [[0, 0, 1, 0, 0, -1, 0]]
    b_lb = [0]
    conc = Concentration(basis_elements=basis_elements,
                         A_eq=A_eq, b_eq=b_eq, A_lb=A_lb, b_lb=b_lb)
    linked = conc._linked_basis
    assert sum(1 for i, num in enumerate(linked) if num == i) == 1

    A_lb = [[1, 0, 0, 1, 0, -1, 0]]
    b_lb = [0]

    conc = Concentration(basis_elements=basis_elements,
                         A_lb=A_lb, b_lb=b_lb)
    linked = conc._linked_basis
    assert sum(1 for i, num in enumerate(linked) if num == i) == 1


def test_two_of_four_linked_basis():
    basis_elements = [["Al", "Mg"], ["X", "O"], ["Ta", "Se"], ["Ta", "O"]]
    A_eq = [[1, 0, 0, 0, -1, 0, 0, 0]]
    b_eq = [0]
    A_lb = [[0, 0, 1, 0, 0, 0, 0, -1]]
    b_lb = [0]
    conc = Concentration(basis_elements=basis_elements,
                         A_eq=A_eq, b_eq=b_eq, A_lb=A_lb, b_lb=b_lb)
    linked = conc._linked_basis
    assert sum(1 for i, num in enumerate(linked) if num == i) == 2





test_full_range()
fixed_composition()
test_conc_range()
test_fix_Ru_composition()
test_formula_unit1()
test_formula_unit2()
test_formula_unit3()
test_formula_unit4()
test_formula_unit5()
test_formula_unit6()
test_three_interlinked_basis()
test_two_of_four_linked_basis()
