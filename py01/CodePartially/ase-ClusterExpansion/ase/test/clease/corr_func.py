"""Unit tests for the corr function class."""
import os
from ase.clease import CEBulk, CorrFunction, Concentration
from ase.clease.corrFunc import equivalent_deco

db_name = "test_corrfunc.db"
conc_args = {"conc_ratio_min_1": [[1, 0]],
             "conc_ratio_max_1": [[0, 1]]}
basis_elements = [["Au", "Cu", "Si"]]
concentration = Concentration(basis_elements=basis_elements)


bc_setting = CEBulk(crystalstructure="fcc", a=4.05, size=[4, 4, 4],
                    concentration=concentration, db_name=db_name,
                    max_cluster_size=3)


def test_trans_matrix():
    """Check that the MIC distance between atoms are correct."""
    atoms = bc_setting.atoms
    tm = bc_setting.trans_matrix
    ref_dist = atoms.get_distance(0, 1, mic=True)
    for indx in range(len(atoms)):
        dist = atoms.get_distance(indx, tm[indx][1], mic=True)
        assert abs(dist - ref_dist) < 1E-5


def get_mic_dists(atoms, cluster):
    """Get the MIC dist."""
    dists = []
    for indx in cluster:
        dist = atoms.get_distances(indx, cluster, mic=True)
        dists.append(dist)
    return dists


def test_order_indep_ref_indx():
    """Check that the order of the elements are independent of the ref index.

    This does only apply for clusters with only inequivalent
    sites
    """
    for _, clst in bc_setting.cluster_info_given_size(3)[0].items():
        if clst["equiv_sites"]:
            # The cluster contains symmetrically equivalent sites
            # and then this test does not apply
            continue
        cluster = clst["indices"]
        cluster_order = clst["order"]

        init_cluster = [0] + list(cluster[0])
        init_cluster = [init_cluster[indx] for indx in cluster_order[0]]

        # Make sure that when the other indices in init_cluster are ref
        # indices, the order is the same
        for ref_indx in cluster[0]:
            found_cluster = False
            for subcluster, order in zip(cluster, cluster_order):
                new_cluster = [ref_indx]
                for indx in subcluster:
                    trans_indx = bc_setting.trans_matrix[ref_indx][indx]
                    new_cluster.append(trans_indx)

                # Check if all elements are the same
                if sorted(new_cluster) == sorted(init_cluster):
                    new_cluster = [new_cluster[indx] for indx in order]
                    found_cluster = True
                    assert init_cluster == new_cluster
            assert found_cluster


def test_interaction_contribution_symmetric_clusters():
    """Test contribution from symmetric clusters.

    Test that when one atom is surrounded by equal atoms,
    the contribution from all clusters within one category is
    the same.
    """
    from ase.build import bulk
    from ase.clease.tools import wrap_and_sort_by_position

    # Create an atoms object that fits with CEBulk
    atoms = bulk("Au", crystalstructure="fcc", a=4.05)
    atoms = atoms * (6, 6, 6)
    atoms = wrap_and_sort_by_position(atoms)

    # Put an Si atom at the origin
    atoms[0].symbol = "Si"

    # Define decoration numbers and make sure they are different
    deco = [[], [], [0, 1], [0, 1, 1], [0, 1, 1, 1]]
    cf = CorrFunction(bc_setting)
    bf = bc_setting.basis_functions
    for size in range(2, 4):
        clusters = bc_setting.cluster_info_given_size(size)[0]
        for _, clst in clusters.items():
            cluster = clst["indices"]
            orders = clst["order"]
            equiv_sites = clst["equiv_sites"]

            equiv_deco = equivalent_deco(deco[size], equiv_sites)
            if len(equiv_deco) == size:
                # Calculate reference contribution for this cluster category
                indices = [0] + list(cluster[0])
                indices = [indices[indx] for indx in orders[0]]
                ref_sp = 0.0
                counter = 0
                for dec in equiv_deco:
                    ref_sp_temp = 1.0
                    for dec_num, indx in zip(dec, indices):
                        ref_sp_temp *= bf[dec_num][atoms[indx].symbol]
                    counter += 1
                    ref_sp += ref_sp_temp
                ref_sp /= counter

                # Calculate the spin product for this category
                sp, count = \
                    cf._sp_same_shape_deco_for_ref_indx(atoms, 0, clst,
                                                        0, deco[size])
                sp /= count
                assert abs(sp - ref_sp) < 1E-8


def test_supercell_consistency():
    from ase.clease.tools import wrap_and_sort_by_position
    basis_elements=[['Li', 'X'], ['O', 'X']]
    concentration = Concentration(basis_elements=basis_elements)
    db_name_sc = "rocksalt_sc.db"
    setting = CEBulk(crystalstructure='rocksalt',
                     a=4.05,
                     size=[1, 1, 1],
                     concentration=concentration,
                     db_name=db_name_sc,
                     max_cluster_size=3,
                     max_cluster_dia=[7.0, 4.0])
    atoms = setting.atoms.copy()
    cf = CorrFunction(setting)
    cf_dict = cf.get_cf(atoms)

    atoms = wrap_and_sort_by_position(atoms*(4, 3, 2))
    cf_dict_sc = cf.get_cf(atoms)
    for k in cf_dict_sc.keys():
        assert abs(cf_dict[k] - cf_dict_sc[k]) < 1E-6
    os.remove(db_name_sc)

test_trans_matrix()
test_order_indep_ref_indx()
test_interaction_contribution_symmetric_clusters()
test_supercell_consistency()
os.remove(db_name)
