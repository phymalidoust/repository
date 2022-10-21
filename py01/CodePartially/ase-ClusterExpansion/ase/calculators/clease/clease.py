"""Calculator for Cluster Expansion."""
import sys
from copy import deepcopy
import numpy as np
from ase.utils import basestring
from ase.atoms import Atoms
from ase.calculators.calculator import Calculator
from ase.clease import CorrFunction, CEBulk, CECrystal
from ase.clease.corrFunc import equivalent_deco


class MovedIgnoredAtomError(Exception):
    """Raised when ignored atoms is moved."""
    pass


class Clease(Calculator):
    """Class for calculating energy using CLEASE.

    Arguments:
    =========
    setting: CEBulk or BulkSapcegroup object

    cluster_name_eci: dictionary of list of tuples containing
                      cluster names and ECI

    init_cf: (optional) correlation function of init_cf

    logfile: file object or str
        If *logfile* is a string, a file with that name will be opened.
        Use '-' for stdout.
    """

    name = 'CLEASE'
    implemented_properties = ['energy']

    def __init__(self, setting, cluster_name_eci=None, init_cf=None,
                 logfile=None):
        Calculator.__init__(self)

        if not isinstance(setting, (CEBulk, CECrystal)):
            msg = "setting must be CEBulk or CECrystal object."
            raise TypeError(msg)
        self.parameters["eci"] = cluster_name_eci
        self.setting = setting
        self.CF = CorrFunction(setting)
        self.norm_factor = self._generate_normalization_factor()

        # check cluster_name_eci and separate them out
        if isinstance(cluster_name_eci, list) and \
           (all(isinstance(i, tuple) for i in cluster_name_eci) or
                all(isinstance(i, list) for i in cluster_name_eci)):
            self.cluster_names = [tup[0] for tup in cluster_name_eci]
            self.eci = np.array([tup[1] for tup in cluster_name_eci])
        elif isinstance(cluster_name_eci, dict):
            self.cluster_names = []
            self.eci = []
            for cluster_name, eci in cluster_name_eci.items():
                self.cluster_names.append(cluster_name)
                self.eci.append(eci)
            self.eci = np.array(self.eci)
        else:
            msg = "'cluster_name_eci' needs to be either (1) a list of tuples "
            msg += "or (2) a dictionary.\n They can be etrieved by "
            msg += "'get_cluster_name_eci' method in Evaluate class."
            raise TypeError(msg)

        # calculate init_cf or convert init_cf to array
        if init_cf is None:
            self.cf = init_cf
        elif isinstance(init_cf, list):
            if all(isinstance(i, (tuple, list)) for i in init_cf):
                cluster_names = [tup[0] for tup in init_cf]
                # cluster_name_eci and init_cf in the same order
                if cluster_names == self.cluster_names:
                    self.cf = np.array([tup[1] for tup in init_cf],
                                       dtype=float)
                # not in the same order
                else:
                    self.cf = []
                    for name in self.cluster_names:
                        indx = cluster_names.index(name)
                        self.cf.append(init_cf[indx][1])
                    self.cf = np.array(self.cf, dtype=float)
            else:
                self.cf = np.array(init_cf, dtype=float)
        elif isinstance(init_cf, dict):
            self.cf = np.array([init_cf[x] for x in self.cluster_names],
                               dtype=float)
        else:
            raise TypeError("'init_cf' needs to be either (1) a list "
                            "of tuples, (2) a dictionary, or (3) numpy array "
                            "containing correlation function in the same "
                            "order as the 'cluster_name_eci'.")

        if self.cf is not None and len(self.eci) != len(self.cf):
            raise ValueError('length of provided ECIs and correlation '
                             'functions do not match')

        # logfile
        if isinstance(logfile, basestring):
            if logfile == '-':
                logfile = sys.stdout
            else:
                logfile = open(logfile, 'a')
        self.logfile = logfile

        self.energy = None
        # reference atoms for calculating the cf and energy for new atoms
        self.ref_atoms = None
        self.ref_cf = None
        # old atoms for the case where one needs to revert to the previous
        # structure (e.g., Monte Carlo Simulation)
        self.old_atoms = None
        self.old_cf = None
        self.symmetry_group = None
        self.is_backround_index = None
        self.num_si = self._get_num_self_interactions()
        self.dupl_tracker = DuplicationCountTracker(self.setting)

    def set_atoms(self, atoms):
        self.atoms = atoms.copy()
        if self.cf is None:
            self.cf = self.CF.get_cf_by_cluster_names(self.atoms,
                                                      self.cluster_names,
                                                      return_type='array')
        self._copy_current_to_ref()
        self._copy_current_to_old()

        if len(self.setting.atoms) != len(atoms):
            msg = "Passed Atoms object and setting.atoms should have "
            msg += "same number of atoms."
            raise ValueError(msg)
        if not np.allclose(atoms.positions, self.setting.atoms.positions):
            msg = "atomic positions of all atoms in the passed Atoms "
            msg += "object and setting.atoms should be the same. "
            raise ValueError(msg)
        self.symmetry_group = np.zeros(len(atoms), dtype=int)
        for symm, indices in enumerate(self.setting.index_by_trans_symm):
            self.symmetry_group[indices] = symm
        self.is_backround_index = np.zeros(len(atoms), dtype=np.uint8)
        self.is_backround_index[self.setting.background_indices] = 1

    def calculate(self, atoms, properties, system_changes):
        """Calculate the energy of the passed atoms object.

        If accept=True, the most recently used atoms object is used as a
        reference structure to calculate the energy of the passed atoms.
        Returns energy.
        """
        self._check_atoms(atoms)
        Calculator.calculate(self, atoms)
        swapped_indices = self.update_energy()
        self.results['energy'] = self.energy
        self.log()
        if len(swapped_indices) == 0:
            return self.energy

        self._copy_ref_to_old()
        self._copy_current_to_ref()

        return self.energy

    def _copy_current_to_ref(self):
        self.ref_atoms = self.atoms.copy()
        self.ref_cf = deepcopy(self.cf)

    def _copy_current_to_old(self):
        self.old_atoms = self.atoms.copy()
        self.old_cf = deepcopy(self.cf)

    def _copy_ref_to_old(self):
        self.old_atoms = self.ref_atoms.copy()
        self.old_cf = deepcopy(self.ref_cf)

    def restore(self, atoms):
        """Restore the old atoms and correlation functions to the reference."""
        if self.old_atoms is None:
            return
        self.ref_atoms = self.old_atoms.copy()
        self.ref_cf = deepcopy(self.old_cf)
        self.atoms = self.old_atoms.copy()
        self.cf = deepcopy(self.old_cf)
        atoms.numbers = self.old_atoms.numbers

    def update_energy(self):
        """Update correlation function and get new energy."""
        swapped_indices = self.update_cf()
        self.energy = self.eci.dot(self.cf) * len(self.atoms)
        return swapped_indices

    @property
    def indices_of_changed_atoms(self):
        """Return the indices of atoms that have been changed."""
        o_numbers = self.ref_atoms.numbers
        n_numbers = self.atoms.numbers
        check = (n_numbers == o_numbers)
        changed = np.argwhere(check == 0)[:, 0]
        changed = np.unique(changed)
        for index in changed:
            if self.is_backround_index[index] and \
                    self.setting.ignore_background_atoms:
                raise MovedIgnoredAtomError("Atom with index {} is a "
                                            "background atom."
                                            "".format(index))

        return np.unique(changed)

    def get_cf_dict(self):
        """Return the correlation functions as a dict"""
        return dict(self.get_cf_list_tup())

    def get_cf_list_tup(self):
        """Return the correlation function as a list of tuples"""
        return zip(self.cluster_names, self.cf)

    def _symbol_by_index(self, indx):
        return [self.ref_atoms[indx].symbol, self.atoms[indx].symbol]

    def _generate_normalization_factor(self):
        """Return a dictionary with all the normalization factors."""
        norm_fact = {}
        for symm, item in enumerate(self.setting.cluster_info):
            num_atoms = len(self.setting.index_by_trans_symm[symm])
            for name, info in item.items():
                if name not in norm_fact.keys():
                    norm_fact[name] = len(info["indices"]) * num_atoms
                else:
                    norm_fact[name] += len(info["indices"]) * num_atoms
        return norm_fact

    def update_cf(self):
        """Update correlation function based on the reference value."""
        swapped_indices = self.indices_of_changed_atoms
        self.cf = deepcopy(self.ref_cf)
        new_symbs = {}
        # Reset the atoms object
        for indx in swapped_indices:
            new_symbs[indx] = self.atoms[indx].symbol
            self.atoms[indx].symbol = self.ref_atoms[indx].symbol

        for indx in swapped_indices:
            # Swap one index at the time
            self.atoms[indx].symbol = new_symbs[indx]
            for i, name in enumerate(self.cluster_names):
                # find c{num} in cluster type
                n = int(name[1])
                if n == 0:
                    continue

                prefix = (name.rpartition('_')[0])
                dec_str = name.rpartition('_')[-1]
                dec = [int(x) for x in dec_str]

                if n == 1:
                    self.cf[i] = self.CF.get_c1(self.atoms, int(dec_str))
                    continue

                symm = self.symmetry_group[indx]

                if prefix not in self.setting.cluster_info[symm].keys():
                    continue
                cluster = self.setting.cluster_info[symm][prefix]
                count = self.norm_factor[prefix]
                cf_tot = self.cf[i] * count
                cf_change = \
                    self._cf_change_by_indx(indx, cluster, dec) /\
                    self.num_si[prefix]
                self.cf[i] = (cf_tot + (n * cf_change)) / count
        return swapped_indices

    def _cf_change_by_indx(self, ref_indx, cluster, deco):
        """Calculate the change in correlation function based on atomic index.

        This method tracks changes in correaltion function due to change in
        element type for atom with index = ref_indx. Passed trans_list refers
        to the indices of atoms that constitute the cluster
        (after translation).
        """
        symbol = self._symbol_by_index(ref_indx)
        b_f = self.setting.basis_functions
        delta_cf = 0.

        eq_sites = cluster["equiv_sites"]
        indx_order = cluster["order"]
        trans_list = self._translate_indx(ref_indx, cluster["indices"])

        equiv_deco = equivalent_deco(deco, eq_sites)
        for dec in equiv_deco:
            for cluster_indx, order, non_trans in \
                    zip(trans_list, indx_order, cluster["indices"]):
                # NOTE: Here cluster_indx is already translated!
                indices = [ref_indx] + cluster_indx
                indices = [indices[indx] for indx in order]
                cf_new = 1.0
                cf_ref = 1.0
                counter = 0
                for j, indx in enumerate(indices):
                    if indx == ref_indx:
                        counter += 1
                        cf_new *= b_f[dec[j]][symbol[1]]
                        cf_ref *= b_f[dec[j]][symbol[0]]
                    else:
                        cf_new *= b_f[dec[j]][self.atoms[indx].symbol]
                        cf_ref *= b_f[dec[j]][self.atoms[indx].symbol]

                # If self interactions we need to down scale the contribution
                # because this cluster is present fewer times.
                # Example: If indices = [0, 23, 42], this particular clusters
                # exists when 0, 23 and 42 is reference index (total 3 times)
                # if the cluster is [0, 23, 23], this cluster exists only
                # when 0 and 23 is reference index (only 2 times).
                # Hence, we need to down scale the contribution by a factor
                # 2/3.
                scale = len(np.unique(indices))/float(len(indices))

                # We only inspect the effect of the change at one reference
                # index. However, we know that if we looped over all reference
                # indices one given clusters would appear exactly m times, if
                # m is the multiplicity factor. On the other hand when
                # inspecting only one cluster, that may not be the case.
                # Example: If indices = [[0, 1, 1], [0, 1, 1], [1, 0, 0]]
                # the situation is as follows. When 0 is reference index
                # [0, 1, 1] occures to times and [1, 0, 0] occures one time
                # When 1 is reference index [0, 1, 1] will occure one time
                # and [1, 0, 0] will occure two times. Since we only
                # inspect one of the cases, we need to correct for the
                # fact that for one particular reference index the number
                # of occurences can be distributed unevenly.
                scale *= self.dupl_tracker.factor(cluster, non_trans, order)
                delta_cf += scale*(cf_new - cf_ref)/len(equiv_deco)
        return delta_cf

    def _translate_indx(self, ref_indx, indx_list):
        tlist = [[] for _ in range(len(indx_list))]
        for i in range(len(indx_list)):
            for j in range(len(indx_list[i])):
                tlist[i].append(
                    self.setting.trans_matrix[ref_indx][indx_list[i][j]])
        return tlist

    def _check_atoms(self, atoms):
        """Check to see if the passed atoms argument valid.

        This method checks that:
            - atoms argument is Atoms object,
            - atoms has the same size and atomic positions as
                (1) setting.atoms,
                (2) reference Atoms object.
        """
        if not isinstance(atoms, Atoms):
            raise TypeError('Passed argument is not Atoms object')
        if len(self.ref_atoms) != len(atoms):
            raise ValueError('Passed atoms does not have the same size '
                             'as previous atoms')
        if not np.allclose(self.ref_atoms.positions, atoms.positions):
            raise ValueError('Atomic postions of the passed atoms are '
                             'different from init_atoms')

    def log(self):
        """Write energy to log file."""
        if self.logfile is None:
            return True
        self.logfile.write('{}\n'.format(self.energy))
        self.logfile.flush()

    def _get_num_self_interactions(self):
        num_si = {}
        for item in self.setting.cluster_info:
            for name, info in item.items():
                if info["size"] <= 1:
                    num_si[name] = 1.
                    continue
                ref_indx = info["ref_indx"]
                num_int = sum(sub.count(ref_indx) for sub in info["indices"])
                num_equiv_indices = 0
                for sub in info["indices"]:
                    bin_count = np.bincount(sub)
                    bin_count[bin_count > 0] -= 1
                    num_equiv_indices += np.sum(bin_count)
                num_int += num_equiv_indices + len(info["indices"])
                num_si[name] = float(num_int)/len(info["indices"])
                num_si[name] = 1.0
        return num_si


def list2str(array):
    return "-".join(str(x) for x in array)


class DuplicationCountTracker(object):
    """Tracks duplication counts and normalization factors.

    Arguments
    ==========
    cluster_info: list of dicts
        The entire info entry in settings
    """
    def __init__(self, setting):
        self.symm_group = np.zeros(len(setting.atoms), dtype=np.uint8)
        for num, group in enumerate(setting.index_by_trans_symm):
            self.symm_group[group] = num

        self.cluster_info = setting.cluster_info
        self.trans_matrix = setting.trans_matrix
        self.occ_count_all = self._occurence_count_all_symm_groups()
        self._norm_factors = self._get_norm_factors()

    def factor(self, cluster, indices, order):
        """Get the normalization factor to correct for self interactions.

        Arguments
        ===========
        cluster: dict
            Dictionary holding information about the cluster

        indices: list
            Indices of the particular sub cluster

        order: list
            Order of the indices in the sub cluster
        """
        key = self.index_key(cluster["ref_indx"], indices, order,
                             cluster["equiv_sites"])
        return self._norm_factors[cluster["symm_group"]][cluster["name"]][key]

    def index_key(self, ref_index, indices, order, equiv_sites):
        """Return a string representing the key for a given order.

        Arguments
        ==========
        ref_index: int
            Reference index

        indices: list
            List representing the indices in a sub-cluster

        order: list
            Order of the indices in the sub cluster
        """
        index_with_ref = [ref_index] + indices
        srt_indices = [index_with_ref[i] for i in order]
        return list2str(self._order_equiv_sites(equiv_sites, srt_indices))

    def _get_norm_factors(self):
        """Calculate all normalization factors."""
        norm_factors = []
        for item in self.cluster_info:
            factors = self._get_norm_factors_per_symm_group(item)
            norm_factors.append(factors)
        return norm_factors

    def _get_norm_factors_per_symm_group(self, clusters):
        """Get normalization factors per symmetry group.

        Arguments:
        =========
        clusters: dict
            Information dict about all clusters in a symmetry group
        """
        norm_factor = {}
        for name, info in clusters.items():
            occ_count = self.occ_count_all[info["symm_group"]][name]
            norm_factor[name] = self._norm_factor(occ_count, info)
        return norm_factor

    def _occurence_count_all_symm_groups(self):
        occ_count_all = []
        for item in self.cluster_info:
            occ_count = {}
            for name, info in item.items():
                occ_count[name] = self._occurence_count(info)
            occ_count_all.append(occ_count)
        return occ_count_all

    def _occurence_count(self, cluster):
        """Count the number of occurences of each sub-cluster in the cluster

        Arguments:
        =========
        cluster: dict
            A dictionary with info about a particular cluster
        """
        occ_count = {}
        for indices, order in zip(cluster["indices"], cluster["order"]):
            key = self.index_key(cluster["ref_indx"], indices, order,
                                 cluster["equiv_sites"])
            occ_count[key] = occ_count.get(key, 0) + 1
        return occ_count

    def _norm_factor(self, occ_count, cluster):
        norm_count = {}
        for k, v in occ_count.items():
            tot_num = self._total_number_of_occurences(k, cluster["name"])
            num_unique = len(set(k.split("-")))
            norm_count[k] = float(tot_num)/(num_unique*v)
        return norm_count

    def _corresponding_subcluster(self, new_ref_indx, target_cluster, name):
        """Find the corresponding cluster when another index is ref_index."""
        cluster = self.cluster_info[self.symm_group[new_ref_indx]][name]
        for sub, order in zip(cluster["indices"], cluster["order"]):
            indices = [self.trans_matrix[new_ref_indx][indx] for indx in sub]
            indices = [new_ref_indx] + indices
            indices = [indices[j] for j in order]
            indices = self._order_equiv_sites(cluster["equiv_sites"], indices)

            if np.allclose(indices, target_cluster):
                subcluster = [cluster["ref_indx"]] + sub
                subcluster = [subcluster[i] for i in order]
                return self._order_equiv_sites(cluster["equiv_sites"],
                                               subcluster)

        raise RuntimeError("There are no matching subcluster. "
                           "This should never happen and is a bug.")

    def _total_number_of_occurences(self, key, name):
        """Get the total number of occurences."""
        indices = list(map(int, key.split("-")))
        tot_num = 0
        for ref_indx in set(indices):
            corr_cluster = \
                self._corresponding_subcluster(ref_indx, indices, name)
            new_key = list2str(corr_cluster)
            tot_num += \
                self.occ_count_all[self.symm_group[ref_indx]][name][new_key]
        return tot_num

    def show(self):
        """Return a string represenatation."""
        print(self._norm_factors)

    def _order_equiv_sites(self, equiv_sites, ordered_indices):
        """After the indices are ordered, adopt a consistent scheme
           within the equivalent sites."""
        for eq_group in equiv_sites:
            equiv_indices = [ordered_indices[i] for i in eq_group]
            equiv_indices.sort()
            for count, i in enumerate(eq_group):
                ordered_indices[i] = equiv_indices[count]
        return ordered_indices
