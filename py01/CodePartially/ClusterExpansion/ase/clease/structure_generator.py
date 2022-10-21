"""Module for generating new structures."""
import os
import math
from random import choice, getrandbits
from copy import deepcopy
import numpy as np
from numpy.linalg import inv, pinv
from ase.db import connect
from ase.clease import CEBulk, CECrystal, CorrFunction
from ase.clease.tools import wrap_and_sort_by_position
from ase.calculators.clease import Clease
from ase.units import kB
import time


class StructureGenerator(object):
    """Base class for generating new strctures."""

    def __init__(self, setting, atoms, struct_per_gen, init_temp=None,
                 final_temp=None, num_temp=5, num_steps_per_temp=10000):
        if not isinstance(setting, (CEBulk, CECrystal)):
            raise TypeError("setting must be CEBulk or CECrystal "
                            "object")

        self.setting = setting
        self.trans_matrix = setting.trans_matrix
        self.cluster_names = self.setting.cluster_names
        self.corrFunc = CorrFunction(setting)
        self.cfm = self._get_full_cf_matrix()
        self.atoms = wrap_and_sort_by_position(atoms.copy())
        if self._is_valid(atoms):
            if len(atoms) != len(setting.atoms):
                raise ValueError("Passed Atoms has a wrong size.")
        else:
            raise ValueError("concentration of the elements in the provided"
                             " atoms is not consistent with the settings.")

        # eci set to 1 to ensure that all correlation functions are included
        # but the energy produced from this should never be used
        eci = {name: 1. for name in self.cluster_names}
        self.calc = Clease(self.setting, cluster_name_eci=eci)
        self.atoms.set_calculator(self.calc)
        self.output_every = 30
        self.init_temp = init_temp
        self.final_temp = final_temp
        self.temp = self.init_temp
        self.num_temp = num_temp
        self.num_steps_per_temp = num_steps_per_temp
        self.alter_composition = True

        # Structures and correlation function associated with
        # the generated structure so far
        self.generated_structure = None
        self.cf_generated_structure = None

    def _is_valid(self, atoms):
        return self.setting.concentration.is_valid(
                self.setting.index_by_basis, atoms)

    def _reset(self):
        pass

    def generate(self):
        """Generate new structures."""
        # Start
        self._reset()
        if self.init_temp is None or self.final_temp is None:
            self.init_temp, self.final_temp = self._determine_temps()

        if self.init_temp <= self.final_temp:
            raise ValueError("Initial temperature must be higher than final"
                             " temperature")
        self._reset()

        temps = np.logspace(math.log10(self.init_temp),
                            math.log10(self.final_temp),
                            self.num_temp)
        now = time.time()
        change_element = self._has_more_than_one_conc()
        for temp in temps:
            self.temp = temp
            num_accepted = 0
            count = 0
            while count < self.num_steps_per_temp:
                count += 1
                if time.time() - now > self.output_every:
                    acc_rate = float(num_accepted)/count
                    print("Temp: {}. {} of {}. Acc. rate: {}"
                          "".format(temp, count, self.num_steps_per_temp,
                                    acc_rate))
                    now = time.time()

                if bool(getrandbits(1)) and self.alter_composition:
                    # Change element Type
                    if change_element:
                        self._change_element_type()
                    else:
                        continue
                else:
                    indx = self._swap_two_atoms()
                    if self.atoms[indx[0]].symbol == \
                            self.atoms[indx[1]].symbol:
                        continue
                self.atoms.get_potential_energy()

                if self._accept():
                    num_accepted += 1
                else:
                    self.calc.restore(self.atoms)

        self._check_consistency()
        cf = self.corrFunc.get_cf(self.atoms, return_type="dict")
        return self.atoms, cf

    def _accept(self):
        raise NotImplementedError('_accept should be implemented in the '
                                  'inherited class.')

    def _estimate_temp_range(self):
        raise NotImplementedError('_estimate_temp_range should be '
                                  'implemented in the inherited class.')

    def _optimal_structure(self):
        raise NotImplementedError('_optimal_structure shoud be implemented '
                                  'in the inherited class.')

    def _determine_temps(self):
        print("Temperature range not given. "
              "Determining the range automatically.")
        self._reset()
        count = 0
        max_count = 100
        now = time.time()
        # To avoid errors, just set the temperature to an arbitrary file
        self.temp = 10000000.0
        while count < max_count:
            if time.time() - now > self.output_every:
                print("Progress ({}%)".format(100*count/max_count))
                now = time.time()

            if bool(getrandbits(1)) and self.alter_composition:
                # Change element Type
                if self._has_more_than_one_conc():
                    self._change_element_type()
                    count += 1
                else:
                    continue
            else:
                indx = self._swap_two_atoms()
                if self.atoms[indx[0]].symbol == \
                        self.atoms[indx[1]].symbol:
                    continue
                count += 1
            self.atoms.get_potential_energy()

            # By calling accept statistics on the correlation
            # function and variance will be collected
            self._accept()
        init_temp, final_temp = self._estimate_temp_range()
        self.temp = init_temp
        print('init_temp= {}, final_temp= {}'.format(init_temp, final_temp))
        return init_temp, final_temp

    def _swap_two_atoms(self):
        """Swap two randomly chosen atoms."""
        indx = np.zeros(2, dtype=int)
        symbol = [None] * 2

        basis_elements = self.setting.basis_elements
        num_basis = self.setting.num_basis

        # pick fist atom and determine its symbol and type
        while True:
            basis = choice(range(num_basis))
            # a basis with only 1 type of element should not be chosen
            if len(basis_elements[basis]) < 2:
                continue
            indx[0] = choice(self.setting.index_by_basis[basis])
            symbol[0] = self.atoms[indx[0]].symbol
            break
        # pick second atom that occupies the same basis.
        while True:
            indx[1] = choice(self.setting.index_by_basis[basis])
            symbol[1] = self.atoms[indx[1]].symbol
            if symbol[1] in basis_elements[basis]:
                break

        # Swap two elements
        self.atoms[indx[0]].symbol = symbol[1]
        self.atoms[indx[1]].symbol = symbol[0]

        return indx

    def _has_more_than_one_conc(self):
        ranges = self.setting.concentration.get_individual_comp_range()
        for r in ranges:
            if abs(r[1] - r[0]) > 0.01:
                return True
        return False

    def _change_element_type(self):
        """Change the type of element for the atom with a given index.

        If index and replacing element types are not specified, they are
        randomly generated.
        """
        basis_elements = self.setting.basis_elements
        num_basis = self.setting.num_basis
        # ------------------------------------------------------
        # Change the type of element for a given index if given.
        # If index not given, pick a random index
        # ------------------------------------------------------
        while True:
            basis = choice(range(num_basis))
            # a basis with only 1 type of element should not be chosen
            if len(basis_elements[basis]) < 2:
                continue

            indx = choice(self.setting.index_by_basis[basis])
            old_symbol = self.atoms[indx].symbol

            # change element type
            new_symbol = choice(basis_elements[basis])
            if new_symbol != old_symbol:
                self.atoms[indx].symbol = new_symbol

                if self.setting.concentration.is_valid(
                     self.setting.index_by_basis, self.atoms):
                    break
                self.atoms[indx].symbol = old_symbol

    def _check_consistency(self):
        # Check to see if the cf is indeed preserved
        final_cf = \
            self.corrFunc.get_cf_by_cluster_names(self.generated_structure,
                                                  self.calc.cluster_names,
                                                  return_type='dict')
        for k in final_cf:
            if abs(final_cf[k] - self.cf_generated_structure[k]) > 1E-6:
                msg = 'Correlation function changed after simulated annealing'
                raise ValueError(msg)

    def _get_full_cf_matrix(self):
        """Get correlation function of every entry in DB."""
        cfm = []
        db = connect(self.setting.db_name)
        for row in db.select(struct_type='initial'):
            cfm.append([row[x] for x in self.cluster_names])
        cfm = np.array(cfm, dtype=float)
        return cfm


class ProbeStructure(StructureGenerator):
    """Generate probe structures.

    Based on simulated annealing according to the recipe in
    PRB 80, 165122 (2009).

    Arguments:
    =========
    setting: CEBulk or BulkSapcegroup object

    atoms: Atoms object
        initial structure to start the simulated annealing

    struct_per_gen: int
        number of structures to be generated per generation

    init_temp: int or float
        initial temperature (does not represent *physical* temperature)

    final_temp: int or float
        final temperature (does not represent *physical* temperature)

    num_temp: int
        number of temperatures to be used in simulated annealing

    num_steps_per_temp: int
        number of steps per temperature in simulated annealing

    approx_mean_var: bool
        whether or not to use a spherical and isotropical distribution
        approximation scheme for determining the mean variance.
        -'True': Assume a spherical and isotropical distribution of
                 structures in the configurational space.
                 Corresponds to eq.4 in PRB 80, 165122 (2009)
        -'False': Use sigma and mu of eq.3 in PRB 80, 165122 (2009)
                  to characterize the distribution of structures in
                  population.
                  Requires pre-sampling of random structures before
                  generating probe structures.
                  Reads sigma and mu from 'probe_structure-sigma_mu.npz' file.
    """

    def __init__(self, setting, atoms, struct_per_gen, init_temp=None,
                 final_temp=None, num_temp=5, num_steps_per_temp=1000,
                 approx_mean_var=False):

        StructureGenerator.__init__(self, setting, atoms, struct_per_gen,
                                    init_temp, final_temp, num_temp,
                                    num_steps_per_temp)
        self.o_cf = self.calc.cf
        self.o_cfm = np.vstack((self.cfm, self.o_cf))
        self.approx_mean_var = approx_mean_var
        fname = 'probe_structure-sigma_mu.npz'
        if not approx_mean_var:
            if os.path.isfile(fname):
                data = np.load(fname)
                self.sigma = data['sigma']
                self.mu = data['mu']
                self.o_mv = mean_variance(self.o_cfm, self.sigma, self.mu)
            else:
                raise IOError("'{}' not found.".format(fname))
        else:
            self.o_mv = mean_variance_approx(self.o_cfm)
        self.avg_mv = 0.0
        self.num_steps = 0
        self.avg_diff = 0.0

        self.min_mv = None

    def _accept(self):
        """Accept the last change."""
        cfm = np.vstack((self.cfm, self.calc.cf))
        if self.approx_mean_var:
            n_mv = mean_variance_approx(cfm)
        else:
            n_mv = mean_variance(cfm, self.sigma, self.mu)

        # Always accept the first move
        if self.generated_structure is None:
            self.generated_structure = self.atoms.copy()
            self.min_mv = n_mv
            self.cf_generated_structure = deepcopy(self.calc.get_cf_dict())
            return True

        if n_mv < self.o_mv:
            self.generated_structure = self.atoms.copy()
            self.min_mv = n_mv
            self.cf_generated_structure = deepcopy(self.calc.get_cf_dict())

        if self.o_mv > n_mv:
            accept = True
        else:
            accept = np.exp((self.o_mv-n_mv)/self.temp) > np.random.uniform()

        self.avg_diff += abs(n_mv - self.o_mv)
        if accept:
            self.o_mv = n_mv

        self.avg_mv += self.o_mv
        self.num_steps += 1
        return accept

    def _estimate_temp_range(self):
        if self.num_steps == 0:
            return 100000.0, 1.0
        avg_diff = self.avg_diff/self.num_steps
        init_temp = 10 * avg_diff
        final_temp = 0.01 * avg_diff
        return init_temp, final_temp

    @property
    def avg_mean_variance(self):
        if self.num_steps == 0:
            return 0.0
        return self.avg_mv/self.num_steps

    def _reset(self):
        self.avg_mv = 0.0
        self.avg_diff = 0.0
        self.num_steps = 0


class EminStructure(StructureGenerator):
    """Generate minimum energy structure.

    Arguments:
    =========
    setting: CEBulk or BulkSapcegroup object

    atoms: Atoms object
        initial structure to start the simulated annealing

    struct_per_gen: int
        number of structures to be generated per generation

    init_temp: int or float
        initial temperature (does not represent *physical* temperature)

    final_temp: int or float
        final temperature (does not represent *physical* temperature)

    num_temp: int
        number of temperatures to be used in simulated annealing

    num_steps_per_temp: int
        number of steps per temperature in simulated annealing

    cluster_name_eci: dict of list of tuples containing
                      cluster names and ECI
    """
    def __init__(self, setting, atoms, struct_per_gen, init_temp=2000,
                 final_temp=10, num_temp=10, num_steps_per_temp=100000,
                 cluster_name_eci=None):
        StructureGenerator.__init__(self, setting, atoms, struct_per_gen,
                                    init_temp, final_temp, num_temp,
                                    num_steps_per_temp)
        self.alter_composition = False
        self.calc = Clease(self.setting, cluster_name_eci=cluster_name_eci)
        self.atoms.set_calculator(self.calc)
        self.old_energy = None
        self.min_energy = None

    def _accept(self):
        """Accept the last change."""
        new_energy = self.atoms.get_potential_energy()

        # Always accept the first move
        if self.generated_structure is None:
            self.old_energy = new_energy
            self.min_energy = new_energy
            self.generated_structure = self.atoms.copy()
            self.cf_generated_structure = deepcopy(self.calc.get_cf_dict())
            return True

        if new_energy < self.min_energy:
            self.min_energy = new_energy
            self.generated_structure = self.atoms.copy()
            self.cf_generated_structure = deepcopy(self.calc.get_cf_dict())

        if new_energy < self.old_energy:
            self.old_energy = new_energy
            return True

        diff = new_energy - self.old_energy
        kT = kB*self.temp
        accept_move = np.random.uniform() < np.exp(-diff/kT)

        if accept_move:
            self.old_energy = new_energy
        return accept_move


def mean_variance_full(cfm):
    prec = precision_matrix(cfm)
    mv = 0.
    for x in range(cfm.shape[0]):
        mv += cfm[x].dot(prec).dot(cfm[x].T)
    mv = mv / cfm.shape[0]
    return mv


def mean_variance(cfm, sigma, mu):
    prec = precision_matrix(cfm)
    return np.trace(prec.dot(sigma)) + mu.dot(prec).dot(mu.T)


def mean_variance_approx(cfm):
    prec = precision_matrix(cfm)
    return np.trace(prec)


def precision_matrix(cfm):
    try:
        prec = inv(cfm.T.dot(cfm))
    # if inverting matrix leads to a singular matrix, reduce the matrix
    except np.linalg.linalg.LinAlgError:
        prec = pinv(cfm.T.dot(cfm))
    return prec
