"""Definition of ClusterExpansionSetting Class.

This module defines the base-class for storing the settings for performing
Cluster Expansion in different conditions.
"""
import os
from itertools import combinations, product
from copy import deepcopy
import numpy as np
from scipy.spatial import cKDTree as KDTree
from ase.db import connect

from ase.clease.floating_point_classification import FloatingPointClassifier
from ase.clease.tools import (wrap_and_sort_by_position, index_by_position,
                              flatten, sort_by_internal_distances,
                              dec_string, get_unique_name,
                              nested_array2list, get_all_internal_distances,
                              distance_string)
from ase.clease.basis_function import BasisFunction
from ase.clease.template_atoms import TemplateAtoms
from ase.clease.concentration import Concentration


class ClusterExpansionSetting(object):
    """Base class for all Cluster Expansion settings."""

    def __init__(self, size=None, supercell_factor=None, dist_num_dec=3,
                 concentration=None, db_name=None, max_cluster_size=4,
                 max_cluster_dia=None, basis_function='sanchez',
                 skew_threshold=4, ignore_background_atoms=False):
        self.kwargs = {'size': size,
                       'supercell_factor': supercell_factor,
                       'db_name': db_name,
                       'max_cluster_size': max_cluster_size,
                       'max_cluster_dia': deepcopy(max_cluster_dia),
                       'dist_num_dec': dist_num_dec,
                       'ignore_background_atoms': ignore_background_atoms}

        if isinstance(concentration, Concentration):
            self.concentration = concentration
        elif isinstance(concentration, dict):
            self.concentration = Concentration.from_dict(concentration)
        else:
            raise TypeError("concentration has to be either dict or "
                            "instance of Concentration")

        self.kwargs["concentration"] = self.concentration.to_dict()
        self.basis_elements = deepcopy(self.concentration.basis_elements)
        self.num_basis = len(self.basis_elements)
        self.db_name = db_name
        self.dist_num_dec = dist_num_dec
        self.size = size
        self.unit_cell = self._get_unit_cell()
        self._tag_unit_cell()
        self.unit_cell_id = self._store_unit_cell()
        self.float_max_dia, self.float_ang, self.float_dist = \
            self._init_floating_point_classifiers()

        if (supercell_factor is not None and size is None and
                max_cluster_dia is None):
            raise ValueError("'max_cluster_dia' must be specified when "
                             "'supercell_factor' is used. ")

        self.template_atoms = TemplateAtoms(supercell_factor=supercell_factor,
                                            size=size,
                                            skew_threshold=skew_threshold,
                                            unit_cell_id=self.unit_cell_id,
                                            db_name=self.db_name)

        self.max_cluster_size = max_cluster_size
        self.max_cluster_dia = max_cluster_dia
        self.all_elements = sorted([item for row in self.basis_elements for
                                    item in row])
        self.ignore_background_atoms = ignore_background_atoms
        self.background_indices = None
        self.num_elements = len(self.all_elements)
        self.unique_elements = sorted(list(set(deepcopy(self.all_elements))))
        self.num_unique_elements = len(self.unique_elements)
        self.index_by_basis = None

        self.cluster_info = []
        self.index_by_trans_symm = []
        self.ref_index_trans_symm = []
        self.kd_trees = None
        self.template_atoms_uid = 0
        for uid in range(self.template_atoms.num_templates):
            self._set_active_template_by_uid(uid)
        # Set the initial template atoms to 0, which is the smallest cell
        self._set_active_template_by_uid(0)
        self._check_cluster_info_consistency()

        unique_element_no_bkg = self.unique_element_without_background()
        if isinstance(basis_function, BasisFunction):
            if basis_function.unique_elements != unique_element_no_bkg:
                raise ValueError("Unique elements in BasiFunction instance "
                                 "is different from the one in settings")
            self.bf_scheme = basis_function

        elif isinstance(basis_function, str):
            if basis_function.lower() == 'sanchez':
                from ase.clease.basis_function import Sanchez
                self.bf_scheme = Sanchez(unique_element_no_bkg)
            elif basis_function.lower() == 'vandewalle':
                from ase.clease.basis_function import VandeWalle
                self.bf_scheme = VandeWalle(unique_element_no_bkg)
            elif basis_function.lower() == "sluiter":
                from ase.clease.basis_function import Sluiter
                self.bf_scheme = Sluiter(unique_element_no_bkg)
            else:
                msg = "basis function scheme {} ".format(basis_function)
                msg += "is not supported."
                raise ValueError(msg)
        else:
            raise ValueError("basis_function has to be instance of "
                             "BasisFunction or a string")

        self.spin_dict = self.bf_scheme.spin_dict
        self.basis_functions = self.bf_scheme.basis_functions

        if len(self.basis_elements) != self.num_basis:
            raise ValueError("list of elements is needed for each basis")

        if not os.path.exists(db_name):
            self._store_data()
        else:
            self._read_data()

    def unique_element_without_background(self):
        """Remove backgound elements."""
        if not self.ignore_background_atoms:
            return self.unique_elements
        symbs = []
        for indx in self.background_indices:
            symbs.append(self.atoms[indx].symbol)
        symbs = list(set(symbs))

        for atom in self.atoms:
            if atom.index in self.background_indices:
                continue
            if atom.symbol in symbs:
                symbs.remove(atom.symbol)
        filtered = [e for e in self.unique_elements if e not in symbs]
        return filtered

    def _group_index_by_basis(self):
        raise NotImplementedError("Has to be implemented in derived classes!")

    def _store_floating_point_classifiers(self):
        """Store classifiers in a separate DB entry."""
        from ase.atoms import Atoms
        db = connect(self.db_name)
        if sum(1 for row in db.select(name="float_classification")) >= 1:
            # Entry already exists
            return

        placeholder = Atoms()
        data = {
            "max_cluster_dia": self.float_max_dia.toJSON(),
            "angles": self.float_ang.toJSON(),
            "float_dist": self.float_dist.toJSON()
        }
        db.write(placeholder, data=data, name="float_classification")

    def _init_floating_point_classifiers(self):
        """Initialize the floating point classifiers from the DB if they
           exist, otherwise initialize a new one."""

        db = connect(self.db_name)
        try:
            row = db.get(name="float_classification")
            max_dia = row.data["max_cluster_dia"]
            angles = row.data["angles"]
            dists = row.data["float_dist"]
            float_max_dia = FloatingPointClassifier.fromJSON(max_dia)
            float_ang = FloatingPointClassifier.fromJSON(angles)
            dists = FloatingPointClassifier.fromJSON(dists)
        except KeyError:
            float_max_dia = FloatingPointClassifier(self.dist_num_dec)
            float_ang = FloatingPointClassifier(0)
            dists = FloatingPointClassifier(self.dist_num_dec)
        return float_max_dia, float_ang, dists

    def _size2string(self):
        """Convert the current size into a string."""
        return "x".join((str(item) for item in self.size))

    def _set_active_template_by_uid(self, uid):
        """Set a fixed template atoms object as the active."""
        self.template_atoms_uid = uid
        self.atoms, self.size = \
            self.template_atoms.get_atoms(uid, return_size=True)
        self.atoms = wrap_and_sort_by_position(self.atoms)
        self.unit_cell_id = self.template_atoms.get_unit_cell_id(uid)

        # Load the corresponding unit cell from db
        db = connect(self.db_name)
        self.unit_cell = db.get(id=self.unit_cell_id).toatoms()

        self.index_by_basis = self._group_index_by_basis()
        self.cluster_info = []

        self.max_cluster_dia, self.supercell_scale_factor = \
            self._get_max_cluster_dia_and_scale_factor(self.max_cluster_dia)

        self.background_indices = self._get_background_indices()
        self.index_by_trans_symm = self._group_indices_by_trans_symmetry()
        self.num_trans_symm = len(self.index_by_trans_symm)
        self.ref_index_trans_symm = [i[0] for i in self.index_by_trans_symm]

        # Read information from database
        # Note that if the data is not found, it will generate
        # the nessecary data structures and store them in the database
        self._read_data()

    def set_active_template(self, size=None, atoms=None, unit_cell_id=None,
                            generate_template=False):
        """Set a new template atoms object."""
        if size is not None and atoms is not None:
            raise ValueError("Specify either size or pass Atoms object.")
        if size is not None:
            uid = self.template_atoms.get_uid_with_given_size(
                size=size, unit_cell_id=unit_cell_id,
                generate_template=generate_template)
        elif atoms is not None:
            uid = self.template_atoms.get_uid_matching_atoms(
                atoms=atoms, generate_template=generate_template)
        else:
            uid = self.template_atoms.weighted_random_template()
        self._set_active_template_by_uid(uid)

    def _tag_unit_cell(self):
        """Add a tag to all the atoms in the unit cell to track the index."""
        for atom in self.unit_cell:
            atom.tag = atom.index

    def _store_unit_cell(self):
        """Store unit cell to the database."""
        db = connect(self.db_name)
        shape = self.unit_cell.get_cell_lengths_and_angles()
        for row in db.select(name='unit_cell'):
            uc_shape = row.toatoms().get_cell_lengths_and_angles()
            if np.allclose(shape, uc_shape):
                return row.id

        uid = db.write(self.unit_cell, name='unit_cell')
        return uid

    def _get_unit_cell(self):
        raise NotImplementedError("This function has to be implemented in "
                                  "in derived classes.")

    def _get_max_cluster_dia_and_scale_factor(self, max_cluster_dia):
        cell = self.atoms.get_cell().T
        min_length = self._get_max_cluster_dia(cell)

        # ------------------------------------- #
        # Get max_cluster_dia in an array form #
        # ------------------------------------- #
        # max_cluster_dia is list or array
        if isinstance(max_cluster_dia, (list, np.ndarray)):
            # Length should be either max_cluster_size+1 or max_cluster_size-1
            if len(max_cluster_dia) == self.max_cluster_size + 1:
                for i in range(2):
                    max_cluster_dia[i] = 0.
                max_cluster_dia = np.array(max_cluster_dia, dtype=float)
            elif len(max_cluster_dia) == self.max_cluster_size - 1:
                max_cluster_dia = np.array(max_cluster_dia, dtype=float)
                max_cluster_dia = np.insert(max_cluster_dia, 0, [0., 0.])
            else:
                raise ValueError("Invalid length for max_cluster_dia.")
        # max_cluster_dia is int or float
        elif isinstance(max_cluster_dia, (int, float)):
            max_cluster_dia *= np.ones(self.max_cluster_size - 1, dtype=float)
            max_cluster_dia = np.insert(max_cluster_dia, 0, [0., 0.])
        # Case for *None* or something else
        else:
            max_cluster_dia = np.ones(self.max_cluster_size - 1, dtype=float)
            max_cluster_dia *= min_length
            max_cluster_dia = np.insert(max_cluster_dia, 0, [0., 0.])

        # --------------------------------- #
        # Get scale_factor in an array form #
        # --------------------------------- #
        atoms = self.atoms

        # TODO: Do we need to do something here?
        # It is not the cell vectors that needs to be twice as long as the
        # maximum cluster distance, but the smallest length of the vector
        # formed by any linear combination of the cell vectors
        lengths = atoms.get_cell_lengths_and_angles()[:3] / 2.0
        scale_factor = max(max_cluster_dia) / lengths
        scale_factor = np.ceil(scale_factor).astype(int)
        scale_factor = self._get_scale_factor(cell, max(max_cluster_dia))
        return max_cluster_dia.round(decimals=self.dist_num_dec), scale_factor

    def _get_max_cluster_dia(self, cell, ret_weights=False):
        lengths = []
        weights = []
        for w in product([-1, 0, 1], repeat=3):
            vec = cell.dot(w)
            if w == (0, 0, 0):
                continue
            lengths.append(np.sqrt(vec.dot(vec)))
            weights.append(w)

        # Introduce tolerance to make max distance strictly
        # smaller than half of the shortest cell dimension
        tol = 2 * 10**(-self.dist_num_dec)
        min_length = min(lengths) / 2.0
        min_length = min_length.round(decimals=self.dist_num_dec) - tol

        if ret_weights:
            min_indx = np.argmin(lengths)
            return min_length, weights[min_indx]
        return min_length

    def _get_scale_factor(self, cell, max_cluster_dia):
        """Compute the scale factor nessecary to resolve max_cluster_dia."""
        cell_to_small = True
        scale_factor = [1, 1, 1]
        orig_cell = cell.copy()
        while cell_to_small:

            # Check what the maximum cluster distance is for the current
            # cell
            max_size, w = self._get_max_cluster_dia(cell, ret_weights=True)
            if max_size < max_cluster_dia:
                # Find which vectors formed the shortest diagonal
                indices_in_w = [i for i, weight in enumerate(w) if weight != 0]
                shortest_vec = -1
                shortest_length = 1E10

                # Find the shortest vector of the ones that formed the
                # shortest diagonal
                for indx in indices_in_w:
                    vec = cell[:, indx]
                    length = np.sqrt(vec.dot(vec))
                    if length < shortest_length:
                        shortest_vec = indx
                        shortest_length = length

                # Increase the scale factor in the direction of the shortest
                # vector in the diagonal by 1
                scale_factor[shortest_vec] += 1

                # Update the cell to the new scale factor
                for i in range(3):
                    cell[:, i] = orig_cell[:, i] * scale_factor[i]
            else:
                cell_to_small = False
        return scale_factor

    def _get_background_indices(self):
        """Get indices of the background atoms."""
        # check if any basis consists of only one element type
        basis = [i for i, b in enumerate(self.basis_elements) if len(b) == 1]

        bkg_indices = []
        for b_indx in basis:
            bkg_indices += self.index_by_basis[b_indx]
        return bkg_indices

    def _get_atoms(self):
        """Create atoms with a user-specified size."""
        atoms = self.unit_cell.copy() * self.size
        return wrap_and_sort_by_position(atoms)

    def _create_kdtrees(self, atoms):
        kd_trees = []
        trans = []

        cell = atoms.get_cell().T
        weights = [-1, 0, 1]
        for comb in product(weights, repeat=3):
            vec = cell.dot(comb) / 2.0
            trans.append(vec)

        # NOTE: If the error message
        # 'The correlation function changed after simulated annealing'
        # appears when probestructures are generated, uncommenting
        # the next line might be a quick fix. However, this introduce
        # a lot of overhead. For big systems one might easily run out of
        # memory.
        # trans += [atom.position for atom in atoms]

        for t in trans:
            shifted = atoms.copy()
            shifted.translate(t)
            shifted.wrap()
            kd_trees.append(KDTree(shifted.get_positions()))
        return kd_trees

    def _group_indices_by_trans_symmetry(self):
        """Group indices by translational symmetry."""
        indices = [a.index for a in self.unit_cell]
        ref_indx = indices[0]
        # Group all the indices together if its atomic number and position
        # sequences are same
        indx_by_equiv = []
        shifted = self.unit_cell.copy()
        shifted = wrap_and_sort_by_position(shifted)
        an = shifted.get_atomic_numbers()
        pos = shifted.get_positions()

        temp = [[indices[0]]]
        equiv_group_an = [an]
        equiv_group_pos = [pos]
        for indx in indices[1:]:
            vec = self.unit_cell.get_distance(indx, ref_indx, vector=True)
            shifted = self.unit_cell.copy()
            shifted.translate(vec)
            shifted = wrap_and_sort_by_position(shifted)
            an = shifted.get_atomic_numbers()
            pos = shifted.get_positions()

            for equiv_group in range(len(temp)):
                if (an == equiv_group_an[equiv_group]).all() and\
                        np.allclose(pos, equiv_group_pos[equiv_group]):
                    temp[equiv_group].append(indx)
                    break
                else:
                    if equiv_group == len(temp) - 1:
                        temp.append([indx])
                        equiv_group_an.append(an)
                        equiv_group_pos.append(pos)

        for equiv_group in temp:
            indx_by_equiv.append(equiv_group)

        # Now we have found the translational symmetry group of all the atoms
        # in the unit cell, now put all the indices of self.atoms into the
        # matrix based on the tag
        indx_by_equiv_all_atoms = [[] for _ in range(len(indx_by_equiv))]
        symm_group_of_tag = [-1 for _ in range(len(self.unit_cell))]
        for gr_id, group in enumerate(indx_by_equiv):
            for item in group:
                symm_group_of_tag[item] = gr_id

        for atom in self.atoms:
            symm_gr = symm_group_of_tag[atom.tag]
            indx_by_equiv_all_atoms[symm_gr].append(atom.index)
        return indx_by_equiv_all_atoms

    def _assign_correct_family_identifier(self):
        """Make the familily IDs increase size."""
        new_names = {}
        prev_id = {}
        for name in self.cluster_family_names_by_size:
            if name == "c0" or name == "c1":
                new_names[name] = name
            else:
                prefix = name.rpartition("_")[0]
                new_id = prev_id.get(prefix, -1) + 1
                new_name = prefix + "_{}".format(new_id)
                new_names[name] = new_name
                prev_id[prefix] = new_id

        new_cluster_info = []
        for item in self.cluster_info:
            new_dict = {}
            for name, info in item.items():
                new_dict[new_names[name]] = info
                new_dict[new_names[name]]["name"] = new_names[name]
            new_cluster_info.append(new_dict)
        self.cluster_info = new_cluster_info

    def _corresponding_indices(self, indices, supercell):
        """Find the indices in supercell that correspond to the ones in
           self.atoms

        Arguments
        ===========
        indices: list of int
            Indices in self.atoms

        supercell: Atoms
            Supercell object where we want to find the indices
            corresponding to the position in self.atoms
        """
        supercell_indices = []
        for indx in indices:
            pos = self.atoms[indx].position
            dist = supercell.get_positions() - pos
            lengths_sq = np.sum(dist**2, axis=1)
            supercell_indices.append(np.argmin(lengths_sq))
        return supercell_indices

    def _create_cluster_information(self):
        """Create a set of parameters describing the structure.

        Return cluster_info.

        cluster_info: list
            list of dictionaries with information of all clusters
            The dictionaries have the following form:
            {
             "name": Unique name describing the cluster.
                     Example:
                        "c3_3p725_1"
                     means it is a 3-body cluster (c3) with a cluster diameter
                     3.725 angstroms (3p275). The last number is a unique
                     family identification number assigned to all cluster
                     families.

             "descriptor": A string that contains a description of the cluster
                           including all of the internal distances and angles.

             "size": Number of atoms in the clusters.

             "symm_group": Translational symmetry group of the cluster

             "ref_indx": Index of a reference atom for the prototype cluster.

             "indices": List containing the indices of atoms in a cluster.
                        There can be more than on set of indices that form the
                        same cluster family, so it is in a list of list format.
                        An example of a three body cluster:
                            ref_indx: 0
                            indices: [[1, 2, 6], [7, 8, 27], [10, 19, 30]]
                        A full list of indices in the cluster is obtained by
                            [ref_indx] + [10, 19, 30] --> [0, 10, 19, 30]

             "order": The order in which the atoms in the clusters are
                      represented. The indices of atoms in a cluster need to be
                      sorted in a prefined order to ensure a consistent
                      assignment of the basis function.
                      With the same 4-body cluster above, the 3 sets of indices
                      can have the order defined as:
                        [[0, 1, 2, 3], [1, 3, 2, 0], [2, 0, 3, 1]].
                      Then, the third cluster in the example is sorted as
                        Unordered: [ref_indx] + [10, 19, 30] -> [0, 10, 19, 30]
                        Ordered: [19, 0, 30, 10]

             "equiv_sites": List of indices of symmetrically equivalent sites.
                            After ordering, the symmetrically equivalent sites
                            are the same for all clusters in the family.
                            The same 4-body cluster example:
                            1) If the clusters have no equivalent sites,
                               this equiv_sites = []
                            2) If atoms in position 1 and 2 of the ordered
                               cluster list are equivalent, then
                               equiv_sites = [[1, 2]]. Which means that
                               0 and 30 in the cluster [19, 0, 30, 10] are
                               equivalent
                            3) If atom 1, 2 are equivalent and atom 0, 3
                               are equivalent equiv_sites = [[1, 2], [0, 3]]
                               For the cluster [19, 0, 30, 10] that means that
                               0 and 30 are equivalent, and 19 and 10 are
                               equivalent
                            4) If all atoms are symmetrically equivalent
                               equiv_sites = [[0, 1, 2, 3]]
            }
        """
        atoms_cpy = self.atoms.copy()
        for atom in atoms_cpy:
            atom.tag = atom.index
        supercell = atoms_cpy*self.supercell_scale_factor
        supercell = wrap_and_sort_by_position(supercell)

        # If the template atoms is not repeated we need to scale it to at
        # least 2x2x2 when all internal distances are extracted
        unit_cell_lengths = self.unit_cell.get_cell_lengths_and_angles()[:3]
        sc_lengths = supercell.get_cell_lengths_and_angles()[:3]
        ratio = np.round(sc_lengths/unit_cell_lengths, decimals=0).astype(int)
        dist_sc_scale_factor = np.array([1, 1, 1])
        dist_sc_scale_factor[ratio == 1] = 2

        supercell.info['distances'] = get_all_internal_distances(
            supercell*dist_sc_scale_factor, max(self.max_cluster_dia))
        kdtrees = self._create_kdtrees(supercell)

        cluster_info = []
        fam_identifier = []

        # Need newer version
        if np.version.version <= '1.13':
            raise ValueError('Numpy version > 1.13 required')

        # determine cluster information for each inequivalent site
        # (based on translation symmetry)
        ref_indices = self._corresponding_indices(
            self.ref_index_trans_symm, supercell)
        # for site, ref_indx in enumerate(self.ref_index_trans_symm):
        for site, ref_indx in enumerate(ref_indices):
            if (supercell[ref_indx].tag in self.background_indices and
                    self.ignore_background_atoms):
                cluster_info.append({})
                continue
            cluster_info_symm = {}
            cluster_info_symm['c0'] = {
                "indices": [],
                "equiv_sites": [],
                "order": [],
                "ref_indx": self.ref_index_trans_symm[site],
                "symm_group": site,
                "descriptor": "empty",
                "name": "c0",
                "max_cluster_dia": 0.0,
                "size": 0
            }

            cluster_info_symm['c1'] = {
                "indices": [],
                "equiv_sites": [],
                "order": [0],
                "ref_indx": self.ref_index_trans_symm[site],
                "symm_group": site,
                "descriptor": "point_cluster",
                "name": 'c1',
                "max_cluster_dia": 0.0,
                "size": 1
            }

            for size in range(2, self.max_cluster_size + 1):
                indices = self.indices_of_nearby_atom(ref_indx, size, kdtrees)
                if self.ignore_background_atoms:
                    indices = [i for i in indices if
                               supercell[i].tag not in self.background_indices]
                indx_set = []
                descriptor_str = []
                order_set = []
                equiv_sites_set = []
                max_cluster_diameter = []
                for k in combinations(indices, size - 1):
                    d = self.get_min_distance((ref_indx,) + k, kdtrees)
                    if max(d) > self.max_cluster_dia[size]:
                        continue
                    order, eq_sites, string_description = \
                        sort_by_internal_distances(supercell, (ref_indx,) + k,
                                                   self.float_ang)
                    descriptor_str.append(string_description)
                    indx_set.append(k)
                    order_set.append(order)
                    equiv_sites_set.append(eq_sites)
                    max_cluster_diameter.append(max(d))

                if not descriptor_str:
                    msg = "There is no cluster with size {}.\n".format(size)
                    msg += "Reduce max_cluster_size or "
                    msg += "increase max_cluster_dia."
                    raise ValueError(msg)

                # categorize the distances
                unique_descriptors = list(set(descriptor_str))
                unique_descriptors = sorted(unique_descriptors, reverse=True)

                for descr in unique_descriptors:
                    if descr not in fam_identifier:
                        fam_identifier.append(descr)

                for desc in unique_descriptors:
                    # Find the maximum cluster diameter of this category
                    indx = descriptor_str.index(desc)
                    max_dia = distance_string(supercell.info["distances"],
                                              max_cluster_diameter[indx])

                    fam_id = fam_identifier.index(desc)
                    name = get_unique_name(size, max_dia, fam_id)

                    cluster_info_symm[name] = {
                        "indices": [],
                        "equiv_sites": equiv_sites_set[indx],
                        "order": [],
                        "ref_indx": self.ref_index_trans_symm[site],
                        "symm_group": site,
                        "descriptor": desc,
                        "name": name,
                        "max_cluster_dia": max_cluster_diameter[indx],
                        "size": size,
                    }

                for x in range(len(indx_set)):
                    category = unique_descriptors.index(descriptor_str[x])
                    max_dia = distance_string(supercell.info["distances"],
                                              max_cluster_diameter[x])
                    fam_id = fam_identifier.index(unique_descriptors[category])
                    name = get_unique_name(size, max_dia, fam_id)

                    sc_index_set = indx_set[x]
                    index_set = []
                    for indx in sc_index_set:
                        index_set.append(int(supercell[indx].tag))
                    cluster_info_symm[name]["indices"].append(index_set)
                    cluster_info_symm[name]["order"].append(order_set[x])

                    assert cluster_info_symm[name]["equiv_sites"] \
                        == equiv_sites_set[x]
                    assert cluster_info_symm[name]["descriptor"] == \
                        descriptor_str[x]

            cluster_info.append(cluster_info_symm)
        self.cluster_info = cluster_info
        self._assign_correct_family_identifier()
        self._store_floating_point_classifiers()

    @property
    def unique_indices(self):
        """Creates a list with the unique indices."""
        all_indices = deepcopy(self.ref_index_trans_symm)
        for item in self.cluster_info:
            for _, info in item.items():
                    all_indices += flatten(info["indices"])
        return list(set(all_indices))

    @property
    def multiplicity_factor(self):
        """Return the multiplicity factor of each cluster."""
        names = self.cluster_family_names
        mult_factor = {name: 0. for name in names}
        name_found = {name: False for name in names}
        normalization = {name: 0 for name in names}
        for name in names:
            for item in self.cluster_info:
                if name not in item.keys():
                    continue
                name_found[name] = True
                cluster = item[name]
                num_in_group = \
                    len(self.index_by_trans_symm[cluster["symm_group"]])
                mult_factor[name] += len(cluster["indices"]) * num_in_group
                normalization[name] += num_in_group

        for name in mult_factor.keys():
            mult_factor[name] = mult_factor[name] / normalization[name]
        for _, found in name_found.items():
            assert found
        return mult_factor

    def cluster_info_by_name(self, name):
        """Get info entries of all clusters with name."""
        name = str(name)
        info = []
        for item in self.cluster_info:
            if name in item.keys():
                info.append(item[name])
        return info

    def _create_translation_matrix(self):
        """Create and return translation matrix.

        Translation matrix maps the indices of the atoms when an atom with the
        translational symmetry is moved to the location of the reference atom.
        Used in conjunction with cluster_indx to calculate the correlation
        function of the structure.
        """
        natoms = len(self.atoms)
        unique_indices = self.unique_indices

        tm = [{} for _ in range(natoms)]

        for i, ref_indx in enumerate(self.ref_index_trans_symm):
            indices = index_by_position(self.atoms)
            tm[ref_indx] = {col: indices[col] for col in unique_indices}

            for indx in self.index_by_trans_symm[i]:
                if indx == ref_indx:
                    continue
                shifted = self.atoms.copy()
                vec = self.atoms.get_distance(indx, ref_indx, vector=True)
                shifted.translate(vec)
                shifted.wrap()

                indices = index_by_position(shifted)
                tm[indx] = {col: indices[col] for col in unique_indices}
        return tm

    def get_min_distance(self, cluster, kd_trees):
        """Get minimum distances.

        Get the minimum distances between the atoms in a cluster according to
        dist_matrix and return the sorted distances (reverse order)
        """
        d = []
        for _, tree in enumerate(kd_trees):
            row = []
            for x in combinations(cluster, 2):
                x0 = tree.data[x[0], :]
                x1 = tree.data[x[1], :]
                row.append(self._get_distance(x0, x1))
            d.append(sorted(row, reverse=True))
        return np.array(min(d))

    def _get_distance(self, x0, x1):
        """Compute the Euclidean distance between two points."""
        diff = x1 - x0
        length = np.sqrt(diff.dot(diff))
        return length

    def indices_of_nearby_atom(self, ref_indx, size, kd_trees):
        """Return the indices of the atoms nearby.

        Indices of the atoms are only included if distances smaller than
        specified by max_cluster_dia from the reference atom index.
        """
        nearby_indices = []
        for tree in kd_trees:
            x0 = tree.data[ref_indx, :]
            result = tree.query_ball_point(x0, self.max_cluster_dia[size])
            nearby_indices += list(result)

        nearby_indices = list(set(nearby_indices))
        nearby_indices.remove(ref_indx)
        return nearby_indices

    @property
    def cluster_family_names(self):
        """Return a list of all cluster names."""
        families = []
        for item in self.cluster_info:
            families += list(item.keys())
        return list(set(families))

    @property
    def cluster_family_names_by_size(self):
        """Return a list of cluster familes sorted by size."""
        sort_list = []
        for item in self.cluster_info:
            for cname, c_info in item.items():
                # Sorted by 1) Number of atoms in cluster
                # 2) diameter (since cname starts with c3_04nn etc.)
                # 3) Unique descriptor
                sort_list.append((c_info["size"],
                                  cname.rpartition("_")[0],
                                  c_info["descriptor"],
                                  cname))
        sort_list.sort()
        sorted_names = []
        for item in sort_list:
            if item[3] not in sorted_names:
                sorted_names.append(item[3])
        return sorted_names

    @property
    def cluster_names(self):
        """Return the cluster names including decoration numbers."""
        names = ["c0"]
        bf_list = list(range(len(self.basis_functions)))
        for item in self.cluster_info:
            for name, info in item.items():
                if info["size"] == 0:
                    continue
                eq_sites = list(info["equiv_sites"])
                for dec in product(bf_list, repeat=info["size"]):
                    dec_str = dec_string(dec, eq_sites)
                    names.append(name + '_' + dec_str)
        return list(set(names))

    def cluster_info_given_size(self, size):
        """Get the cluster info of all clusters with a given size."""
        clusters = []
        for item in self.cluster_info:
            info_dict = {}
            for key, info in item.items():
                if info["size"] == size:
                    info_dict[key] = info
            clusters.append(info_dict)
        return clusters

    def _store_data(self):
        size_str = "x".join(str(s) for s in self.size)
        num = self.template_atoms_uid
        num_templates = self.template_atoms.num_templates
        print('Generating cluster data for template with size: {}. '
              '({} of {})'.format(size_str, num+1, num_templates))

        self._create_cluster_information()
        self.trans_matrix = self._create_translation_matrix()
        db = connect(self.db_name)
        data = {'cluster_info': self.cluster_info,
                'trans_matrix': self.trans_matrix}
        try:
            row = db.get(name="template", size=self._size2string(),
                         unit_cell_id=self.unit_cell_id)
            db.update(row.id, data=data)
        except KeyError:
            db.write(self.atoms, name='template', data=data,
                     size=self._size2string(), unit_cell_id=self.unit_cell_id)

    def _read_data(self):
        db = connect(self.db_name)
        try:
            select_cond = [('name', '=', 'template'),
                           ('size', '=', self._size2string()),
                           ('unit_cell_id', '=', self.unit_cell_id)]
            row = db.get(select_cond)
            self.cluster_info = row.data.cluster_info
            self._info_entries_to_list()
            self.trans_matrix = row.data.trans_matrix
        except KeyError:
            self._store_data()
        except (AssertionError, AttributeError, RuntimeError):
            self.reconfigure_settings()

    def _info_entries_to_list(self):
        """Convert entries in cluster info to list."""
        for info in self.cluster_info:
            for name, cluster in info.items():
                cluster['indices'] = nested_array2list(cluster['indices'])
                cluster['equiv_sites'] = \
                    nested_array2list(cluster['equiv_sites'])
                cluster['order'] = nested_array2list(cluster['order'])

    def _get_name_indx(self, unique_name):
        size = int(unique_name[1])
        for symm in range(self.num_trans_symm):
            name_list = self.cluster_names[symm][size]
            try:
                n_indx = name_list.index(unique_name)
                return symm, n_indx
            except ValueError:
                continue

    def _group_index_by_basis_group(self):
        if self.concentration.grouped_basis is None:
            return self.index_by_basis

        index_by_grouped_basis = []
        for group in self.concentration.grouped_basis:
            indices = []
            for basis in group:
                indices.extend(self.index_by_basis[basis])
            index_by_grouped_basis.append(indices)

        for basis in index_by_grouped_basis:
            basis.sort()
        return index_by_grouped_basis

    def view_clusters(self):
        """Display all clusters along with their names."""
        from ase.gui.gui import GUI
        from ase.gui.images import Images

        # set the active template which has the best chance of displaying
        # the largest diameter
        sizes = []
        for i in range(self.template_atoms.num_templates):
            atoms = self.template_atoms.get_atoms(i)
            lengths = atoms.get_cell_lengths_and_angles()[:3]
            sizes.append(lengths)

        uid = 0
        candidates = []
        for i, size in enumerate(sizes):
            if i == 0:
                candidates.append(min(size))
                uid = i
                continue
            if min(size) <= max(candidates):
                continue
            uid = i
        self._set_active_template_by_uid(uid)
        if max(self.supercell_scale_factor) > 1:
            print("Warning: the largest template atoms in DB is too small to "
                  "accrately display large clusters. \nPlease change the "
                  "'size' or 'supercell_factor' to include the largest "
                  "template.")

        already_included_names = []
        cluster_atoms = []
        for unique_name in self.cluster_family_names_by_size:
            if unique_name in already_included_names:
                continue
            already_included_names.append(unique_name)
            for symm, entry in enumerate(self.cluster_info):
                if unique_name in entry:
                    cluster = entry[unique_name]
                    break
            if cluster["size"] <= 1:
                continue
            ref_indx = self.ref_index_trans_symm[symm]
            name = cluster["name"]

            atoms = self.atoms.copy()

            keep_indx = [ref_indx] + list(cluster["indices"][0])
            equiv = list(cluster["equiv_sites"])
            order = list(cluster["order"][0])

            if order is not None:
                keep_indx = [keep_indx[n] for n in order]

            for tag, indx in enumerate(keep_indx):
                atoms[indx].tag = tag
            if equiv:
                for group in equiv:
                    for i in range(1, len(group)):
                        atoms[keep_indx[group[i]]].tag = \
                            atoms[keep_indx[group[0]]].tag
            # atoms = create_cluster(atoms, keep_indx)

            # Extract the atoms in cluster for visualization
            atoms = atoms[keep_indx]
            positions = atoms.get_positions()
            cell = atoms.get_cell()
            origin = cell[0, :] + cell[1, :] + cell[2, :]
            supercell = atoms.copy()*(3, 3, 3)
            size = len(keep_indx)
            origin += positions[0, :]
            pos_sc = supercell.get_positions()
            lengths_sq = np.sum((pos_sc - origin)**2, axis=1)
            indices = np.argsort(lengths_sq)[:size]
            atoms = supercell[indices]

            atoms.info = {'name': name}
            cluster_atoms.append(atoms)

        images = Images()
        images.initialize(cluster_atoms)
        gui = GUI(images, expr='')
        gui.show_name = True
        gui.run()
        self._set_active_template_by_uid(0)

    def reconfigure_settings(self):
        """Reconfigure templates stored in DB file."""
        db = connect(self.db_name)
        # Reconfigure the float point classification based on setting
        ids = [row.id for row in db.select(name='float_classification')]
        db.delete(ids)
        self.float_max_dia, self.float_ang, self.float_dist = \
            self._init_floating_point_classifiers()

        # Reconfigure the cluster information for each template based on
        # current max_cluster_size and max_cluster_dia
        for uid in range(self.template_atoms.num_templates):
            self._set_active_template_by_uid(uid)
            self._store_data()
        self._set_active_template_by_uid(0)
        print('Cluster data updated for all templates.\n'
              'You should also reconfigure DB entries (in CorrFunction class) '
              'to make the information on each structure to be consistent. ')

    def _check_first_elements(self):
        basis_elements = self.basis_elements
        num_basis = self.num_basis
        # This condition can be relaxed in the future
        first_elements = []
        for elements in basis_elements:
            first_elements.append(elements[0])
        if len(set(first_elements)) != num_basis:
            raise ValueError("First element of different basis should not be "
                             "the same.")

    def save(self, filename):
        """Write Setting object to a file in JSON format.

        Arguments:
        =========
        filename: str
            Name of the file to store the necessary settings to initialize
            the Cluster Expansion calculations.
        """
        class_types = ['CEBulk', 'CECrystal']
        if type(self).__name__ not in class_types:
            raise NotImplementedError('Class {} '
                                      'is not supported.'
                                      ''.format(type(self).__name__))

        import json
        if type(self).__name__ == 'CEBulk':
            self.kwargs['classtype'] = 'CEBulk'
        else:
            self.kwargs['classtype'] = 'CECrystal'
        # Write keyword arguments necessary for initializing the class
        with open(filename, 'w') as outfile:
            json.dump(self.kwargs, outfile, indent=2)

    def _get_shortest_distance_in_unitcell(self):
        """Find the shortest distance between the atoms in the unit cell."""
        if len(self.unit_cell) == 1:
            lengths = self.unit_cell.get_cell_lengths_and_angles()[:3]
            return min(lengths)

        dists = []
        for ref_atom in range(len(self.unit_cell)):
            indices = list(range(len(self.unit_cell)))
            indices.remove(ref_atom)
            dists += list(self.unit_cell.get_distances(ref_atom, indices,
                                                       mic=True))
        return min(dists)

    def _check_cluster_info_consistency(self):
        """Check that cluster names in all templates' info entries match."""
        db = connect(self.db_name)
        names = []
        descriptors = []
        equiv_sites = {}
        mult_factor = {}
        for row in db.select(name='template'):
            cluster_info = row.data["cluster_info"]

            new_names = []
            new_desc = []
            new_equiv_sites = {}
            new_mult_factor = {}
            # Extract all names
            for item in cluster_info:
                for k, v in item.items():
                    new_names.append(k)
                    new_desc.append(v["descriptor"])
                    if len(v["equiv_sites"]) == 0:
                        new_equiv_sites[k] = 0
                    else:
                        new_equiv_sites[k] = len(v["equiv_sites"][0])
                    new_mult_factor[k] = len(v["indices"])

            new_names = sorted(new_names)
            new_desc = sorted(new_desc)

            if not names:
                # This is the first entry
                names = deepcopy(new_names)
                descriptors = deepcopy(new_desc)
                equiv_sites = deepcopy(new_equiv_sites)
                mult_factor = deepcopy(new_mult_factor)

            assert new_names == names
            assert descriptors == new_desc
            assert equiv_sites == new_equiv_sites
            assert mult_factor == new_mult_factor

    def subclusters(self, cluster_name):
        """Return all the sub-clusters of cluster."""
        info = self.cluster_info_by_name(cluster_name)

        sub_clst = []

        # Loop over symmetry groups
        for symm, item in enumerate(info):
            indices = item["indices"]
            size = item["size"]

            # Search in clusters up it the current size
            for s in range(size):
                clst_size = self.cluster_info_given_size(s)[symm]
                for k, v in clst_size.items():
                    if self._is_subcluster(v["indices"], indices):
                        sub_clst.append(v["name"])
        return list(set(sub_clst))

    def _is_subcluster(self, small_cluster, large_cluster):
        """Return True if small cluster is part of large cluster."""
        if len(small_cluster) == 0:
            # Singlet and empty is always a sub cluster
            return True

        if len(small_cluster[0]) >= len(large_cluster[0]):
            raise ValueError("A cluster with size {} cannot be a "
                             "subcluster of a cluster with size {}"
                             "".format(len(small_cluster[0],
                                       len(large_cluster[0]))))

        return any(set(s1).issubset(s2) for s1 in small_cluster
                   for s2 in large_cluster)
        
