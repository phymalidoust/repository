import os
from ase.clease import CEBulk, CECrystal, Concentration

db_name = 'sf.db'

def test_fcc():
    conc = Concentration(basis_elements=[['Au', 'Cu']])

    setting = CEBulk(crystalstructure='fcc', a=4.0, 
                     cubic=True, supercell_factor=8,
                     concentration=conc, db_name=db_name, max_cluster_size=4,
                     max_cluster_dia=4.0, basis_function='sanchez',
                     skew_threshold=4, ignore_background_atoms=False)
    
    assert setting.template_atoms.num_templates == 6
    
    os.remove(db_name)

    setting = CEBulk(crystalstructure='fcc', a=4.0,
                     cubic=True, supercell_factor=8, size=[2, 2, 2],
                     concentration=conc, db_name=db_name, max_cluster_size=4,
                     max_cluster_dia=4.0, basis_function='sanchez',
                     skew_threshold=4, ignore_background_atoms=False)
    
    assert setting.template_atoms.num_templates == 1

    os.remove(db_name)


def test_crystal():
    basis_elements = [['O', 'X'], ['O', 'X'],
                      ['O', 'X'], ['Ta']]
    grouped_basis = [[0, 1, 2], [3]]
    concentration = Concentration(basis_elements=basis_elements,
                                  grouped_basis=grouped_basis)

    setting = CECrystal(basis=[(0., 0., 0.),
                               (0.3894, 0.1405, 0.),
                               (0.201, 0.3461, 0.5),
                               (0.2244, 0.3821, 0.)],
                        spacegroup=55,
                        cellpar=[6.25, 7.4, 3.83, 90, 90, 90],
                        supercell_factor=10,
                        concentration=concentration,
                        db_name=db_name,
                        basis_function='sluiter',
                        max_cluster_size=3,
                        max_cluster_dia=3.0)
    
    assert setting.template_atoms.num_templates == 20

    os.remove(db_name)


test_fcc()
test_crystal()
