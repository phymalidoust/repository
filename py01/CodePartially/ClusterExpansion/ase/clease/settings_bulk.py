"""Definitions of Cluster Expansion settings for bulk.

Cluster Expansion can be peformed on bulk material using either CEBulk
or CECrystal class.
"""
import numpy as np
from ase.build import bulk
from ase.spacegroup import crystal, Spacegroup
from ase.clease.tools import wrap_and_sort_by_position
from ase.clease.settings import ClusterExpansionSetting
from copy import deepcopy


class CEBulk(ClusterExpansionSetting):
    """Store CE settings on bulk materials defined based on crystal structures.

    Arguments:
    =========
    basis_elements: list
        List of chemical symbols of elements to occupy each basis.
        Even for the cases where there is only one basis (e.g., fcc, bcc, sc),
        a list of symbols should be grouped by basis as in [['Cu', 'Au']]
        (note the nested list form).

    crystalstructure: str
        Must be one of sc, fcc, bcc, hcp, diamond, zincblende, rocksalt,
        cesiumchloride, fluorite or wurtzite.

    a: float
        Lattice constant.

    c: float
        Lattice constant.

    covera: float
        c/a ratio used for hcp.  Default is ideal ratio: sqrt(8/3).

    u: float
        Internal coordinate for Wurtzite structure.

    orthorhombic: bool
        Construct orthorhombic unit cell instead of primitive cell
        which is the default.

    cubic: bool
        Construct cubic unit cell if possible.

    size: list
        size of the supercell (e.g., [2, 2, 2] for 2x2x2 cell)

    conc_args: dict
        ratios of the elements for different concentrations.

    db_name: str
        name of the database file

    max_cluster_size: int
        maximum size (number of atoms in a cluster)

    max_cluster_dia: float or int
        maximum diameter of cluster (in angstrom)

    dist_num_dec: int
        number of decimal places used to determine the distances between atoms

    ignore_background_atoms: bool
        if `True`, a basis consisting of only one element type will be ignored
        when creating clusters.
    """

    def __init__(self, crystalstructure=None,
                 a=None, c=None, covera=None, u=None, orthorhombic=False,
                 cubic=False, size=None, supercell_factor=None,
                 concentration=None, db_name=None, max_cluster_size=4,
                 max_cluster_dia=None, basis_function='sanchez',
                 dist_num_dec=3, skew_threshold=4,
                 ignore_background_atoms=False):

        # Initialization
        self.structures = {'sc': 1, 'fcc': 1, 'bcc': 1, 'hcp': 1, 'diamond': 1,
                           'zincblende': 2, 'rocksalt': 2, 'cesiumchloride': 2,
                           'fluorite': 3, 'wurtzite': 2}
        self.crystalstructure = crystalstructure
        self.a = a
        self.c = c
        self.covera = covera
        self.u = u
        self.orthorhombic = orthorhombic
        self.cubic = cubic

        ClusterExpansionSetting.__init__(self, size, supercell_factor,
                                         dist_num_dec, concentration, db_name,
                                         max_cluster_size, max_cluster_dia,
                                         basis_function, skew_threshold,
                                         ignore_background_atoms)

        # Save raw input arguments for save/load. The arguments gets altered
        # during the initalization process to handle 'ignore_background_atoms'
        # case
        self.kwargs.update({'crystalstructure': crystalstructure,
                            'a': a,
                            'c': c,
                            'covera': covera,
                            'u': u,
                            'orthorhombic': orthorhombic,
                            'cubic': cubic})
        num_basis = len(self.concentration.orig_basis_elements)
        if num_basis != self.structures[self.crystalstructure]:
            msg = "{} has {} basis. ".format(
                self.crystalstructure, self.structures[self.crystalstructure])
            msg += "The number of basis specified by basis_elements is "
            msg += "{}".format(num_basis)
            raise ValueError(msg)

        self._check_first_elements()

    def _get_unit_cell(self):
        basis_elements = self.concentration.orig_basis_elements
        num_basis = len(basis_elements)
        if num_basis == 1:
            atoms = bulk(name='{}'.format(basis_elements[0][0]),
                         crystalstructure=self.crystalstructure, a=self.a,
                         c=self.c, covera=self.covera, u=self.u,
                         orthorhombic=self.orthorhombic, cubic=self.cubic)

        elif num_basis == 2:
            atoms = bulk(name='{}{}'.format(basis_elements[0][0],
                                            basis_elements[1][0]),
                         crystalstructure=self.crystalstructure, a=self.a,
                         c=self.c, covera=self.covera, u=self.u,
                         orthorhombic=self.orthorhombic, cubic=self.cubic)

        else:
            atoms = bulk(name='{}{}{}'.format(basis_elements[0][0],
                                              basis_elements[1][0],
                                              basis_elements[2][0]),
                         crystalstructure=self.crystalstructure, a=self.a,
                         c=self.c, covera=self.covera, u=self.u,
                         orthorhombic=self.orthorhombic, cubic=self.cubic)
        atoms = wrap_and_sort_by_position(atoms)
        return atoms

    def _group_index_by_basis(self):
        indx_by_basis = []
        for basis in self.basis_elements:
            indx_by_basis.append([a.index for a in self.atoms if
                                  a.symbol == basis[0]])

        for basis in indx_by_basis:
            basis.sort()
        self.index_by_basis = indx_by_basis
        return self.index_by_basis

    def _group_index_by_basis_group(self):
        return self.index_by_basis

    @staticmethod
    def load(filename):
        """Load settings from a file in JSON format.

        Arguments:
        =========
        filename: str
            Name of the file that has the settings.
        """
        import json
        with open(filename, 'r') as infile:
            kwargs = json.load(infile)
        classtype = kwargs.pop("classtype")
        if classtype != 'CEBulk':
            raise TypeError('Loaded setting file is not for CEBulk class')
        return CEBulk(**kwargs)


class CECrystal(ClusterExpansionSetting):
    """Store CE settings on bulk materials defined based on space group.

    Arguments:
    =========
    basis_elements: list
        List of chemical symbols of elements to occupy each basis.

    basis: list of scaled coordinates
        Positions of the unique sites corresponding to symbols given
        either as scaled positions or through an atoms instance.

    spacegroup: int | string | Spacegroup instance
        Space group given either as its number in International Tables
        or as its Hermann-Mauguin symbol.

    cell: 3x3 matrix
        Unit cell vectors.

    cellpar: [a, b, c, alpha, beta, gamma]
        Cell parameters with angles in degree. Is not used when `cell`
        is given.

    ab_normal: vector
        Is used to define the orientation of the unit cell relative
        to the Cartesian system when `cell` is not given. It is the
        normal vector of the plane spanned by a and b.

    primitive_cell: bool
        Wheter to return the primitive instead of the conventional
        unit cell.

    size: list of 3 positive integers
        How many times the conventional unit cell should be repeated
        in each direction.

    conc_args: dict
        ratios of the elements for different concentrations.

    db_name: str
        name of the database file

    max_cluster_size: int
        maximum size (number of atoms in a cluster)

    max_cluster_dia: float or int
        maximum diameter of cluster (in angstrom)

    dist_num_dec: int
        number of decimal places used to determine the distances between atoms

    ignore_background_atoms: bool
        if `True`, a basis consisting of only one element type will be ignored
        when creating clusters.
    """

    def __init__(self, basis=None, spacegroup=1,
                 cell=None, cellpar=None, ab_normal=(0, 0, 1), size=None,
                 supercell_factor=None, primitive_cell=False,
                 concentration=None, db_name=None, max_cluster_size=4,
                 max_cluster_dia=None, basis_function='sanchez',
                 dist_num_dec=3, skew_threshold=4,
                 ignore_background_atoms=False):

        # Initialization
        self.basis = basis
        self.spacegroup = spacegroup
        self.cell = cell
        self.cellpar = cellpar
        self.ab_normal = ab_normal
        self.primitive_cell = primitive_cell
        self.symbols = []
        num_basis = len(concentration.orig_basis_elements)
        for x in range(num_basis):
            self.symbols.append(concentration.orig_basis_elements[x][0])

        ClusterExpansionSetting.__init__(self, size, supercell_factor,
                                         dist_num_dec, concentration, db_name,
                                         max_cluster_size, max_cluster_dia,
                                         basis_function, skew_threshold,
                                         ignore_background_atoms)

        # Save raw input arguments for save/load. The arguments gets altered
        # during the initalization process to handle 'ignore_background_atoms'
        # case
        self.kwargs.update({'basis': deepcopy(basis),
                            'spacegroup': spacegroup,
                            'cell': cell,
                            'cellpar': cellpar,
                            'ab_normal': ab_normal,
                            'primitive_cell': primitive_cell})

        self._check_first_elements()

    def _get_unit_cell(self):
        atoms = crystal(symbols=self.symbols, basis=self.basis,
                        spacegroup=self.spacegroup, cell=self.cell,
                        cellpar=self.cellpar, ab_normal=self.ab_normal,
                        size=[1, 1, 1], primitive_cell=self.primitive_cell)
        atoms = wrap_and_sort_by_position(atoms)
        return atoms

    def _group_index_by_basis(self):
        num_basis = len(self.concentration.orig_basis_elements)
        indx_by_basis = [[] for _ in range(num_basis)]
        sg = Spacegroup(self.spacegroup)
        sites, kinds = sg.equivalent_sites(self.basis)
        scale_factor = self.size

        # account for the case where a supercell is needed
        if not np.array_equal(scale_factor, np.array([1, 1, 1])):
            sites = np.divide(sites, scale_factor)
            # x dimension
            kinds = np.tile(kinds, scale_factor[0])
            sites_temp = np.copy(sites)
            for x in range(1, scale_factor[0]):
                shift = np.add(sites_temp, [float(x) / scale_factor[0], 0, 0])
                sites = np.append(sites, shift, axis=0)
            # y dimension
            kinds = np.tile(kinds, scale_factor[1])
            sites_temp = np.copy(sites)
            for y in range(1, scale_factor[1]):
                shift = np.add(sites_temp, [0, float(y) / scale_factor[1], 0])
                sites = np.append(sites, shift, axis=0)
            # z dimension
            kinds = np.tile(kinds, scale_factor[2])
            sites_temp = np.copy(sites)
            for z in range(1, scale_factor[2]):
                shift = np.add(sites_temp, [0, 0, float(z) / scale_factor[2]])
                sites = np.append(sites, shift, axis=0)

        positions = self.atoms.get_scaled_positions()
        for i, site in enumerate(sites):
            for j, pos in enumerate(positions):
                # Avoid position to be very close to 1.0 (e.g., 0.99999999)
                pos -= np.isclose(pos, 1.0, atol=1.e-10, rtol=0.0)
                site -= np.isclose(site, 1.0, atol=1.e-10, rtol=0.0)

                indx = None
                if np.allclose(site, pos):
                    indx = j
                    break

            indx_by_basis[kinds[i]].append(indx)

        for basis in indx_by_basis:
            basis.sort()
        self.index_by_basis = indx_by_basis
        return self._group_index_by_basis_group()
        # return indx_by_basis

    @staticmethod
    def load(filename):
        """Load settings from a file in JSON format.

        Arguments:
        =========
        filename: str
            Name of the file that has the settings.
        """
        import json
        with open(filename, 'r') as infile:
            kwargs = json.load(infile)
        classtype = kwargs.pop("classtype")
        if classtype != 'CECrystal':
            raise TypeError('Loaded setting file is not for CEBulk class')
        return CEBulk(**kwargs)
