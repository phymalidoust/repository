"""Test suite for TemplateAtoms."""
import os
from ase.clease.template_atoms import TemplateAtoms
from ase.build import bulk
from ase.db import connect

db_name = 'templates.db'
def test_fcc():
    unit_cell = bulk("Cu")
    db = connect(db_name)
    db.write(unit_cell, name='unit_cell')

    template_atoms = TemplateAtoms(supercell_factor=27, size=None,
                                   skew_threshold=4, unit_cell_id=0,
                                   db_name=db_name)
    dims = template_atoms.get_size()
    ref = [[1, 1, 1], [1, 1, 2], [2, 2, 2], [2, 2, 3], [2, 2, 4],
           [2, 2, 5], [2, 3, 3], [2, 3, 4], [3, 3, 3]]
    assert dims == ref

    os.remove(db_name)


def test_hcp():
    unit_cell = bulk("Mg")
    db = connect(db_name)
    db.write(unit_cell, name='unit_cell')
    template_atoms = TemplateAtoms(supercell_factor=27, size=None,
                                   skew_threshold=5, unit_cell_id=0,
                                   db_name=db_name)
    dims = template_atoms.get_size()
    ref = [[1, 1, 1], [1, 1, 2], [1, 2, 1], [1, 2, 2], [1, 3, 1], [1, 3, 2],
           [1, 4, 1], [2, 2, 1], [2, 2, 2], [2, 2, 3], [2, 2, 4], [2, 2, 5],
           [2, 3, 1], [2, 3, 2], [2, 3, 3], [2, 3, 4], [2, 4, 1], [2, 4, 2],
           [2, 4, 3], [2, 5, 1], [2, 5, 2], [2, 6, 1], [2, 6, 2], [3, 3, 1],
           [3, 3, 2], [3, 3, 3], [3, 4, 1], [3, 4, 2], [3, 5, 1], [3, 6, 1],
           [4, 4, 1], [4, 5, 1]]

    assert dims == ref

    os.remove(db_name)


test_fcc()
test_hcp()
