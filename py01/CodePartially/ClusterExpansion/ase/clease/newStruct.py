"""Module for generating new structures for training."""
import os
import numpy as np
from random import shuffle
from copy import deepcopy
from functools import reduce

from ase.io import read
from ase.db import connect
from ase.atoms import Atoms
from ase.utils.structure_comparator import SymmetryEquivalenceCheck

from ase.clease import CEBulk, CECrystal, CorrFunction
from ase.clease.structure_generator import ProbeStructure, EminStructure
from ase.clease.tools import wrap_and_sort_by_position

try:
    from math import gcd
except ImportError:
    from fractions import gcd

max_attempt = 10


class MaxAttemptReachedError(Exception):
    """Raised when number of try reaches 10."""

    pass


# class GenerateStructures(object):
class NewStructures(object):
    """Generate new structure in Atoms object format.

    Arguments:
    =========
    setting: CEBulk or BulkSapcegroup object

    generation_number: int
        a generation number to be assigned to the newly generated structure.

    struct_per_gen: int
        number of structures to generate per generation.
    """

    def __init__(self, setting, generation_number=None, struct_per_gen=5):
        if not isinstance(setting, (CEBulk, CECrystal)):
            msg = "setting must be CEBulk or CECrystal object"
            raise TypeError(msg)
        self.setting = setting
        self.db = connect(setting.db_name)
        self.corrfunc = CorrFunction(self.setting)
        self.struct_per_gen = struct_per_gen
        if generation_number is None:
            self.gen, self.num_in_gen = self._determine_gen_number()
        else:
            self.gen = generation_number
            self.num_in_gen = len([row.id for row in
                                   self.db.select(gen=self.gen)])
        if self.num_in_gen >= self.struct_per_gen:
            print("There are {} structures in generation {} in DB and "
                  "struct_per_gen = {}. No more structures generated."
                  .format(self.num_in_gen, self.gen, self.struct_per_gen))
            exit(0)
        self.num_to_gen = self.struct_per_gen - self.num_in_gen

    def generate_probe_structure(self, atoms=None, size=None,
                                 unit_cell_id=None, init_temp=None,
                                 final_temp=None, num_temp=5,
                                 num_steps_per_temp=1000,
                                 approx_mean_var=True,
                                 num_samples_var=10000):
        """Generate a probe structure according to PRB 80, 165122 (2009).

        Arguments:
        =========
        atoms: Atoms object
            Atoms object with the desired cell size and shape of the new
            structure.

        size: list of length 3
            (ignored if atoms is given)
            If specified, the structure with the provided size is generated.
            If None, the size will be generated randomly with a bias towards
                more cubic cells (i.e., cell with similar magnitudes of vectors
                a, b and c)

        unit_cell_id: int
            (only used when size is used)
            The ID of the unit cell in the database to be used

        init_temp: int or float
            initial temperature (does not represent *physical* temperature)

        final_temp: int or float
            final temperature (does not represent *physical* temperature)

        num_temp: int
            number of temperatures to be used in simulated annealing

        num_steps_per_temp: int
            number of steps in simulated annealing

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
                      sigma and mu are generated and stored in
                      'probe_structure-sigma_mu.npz' file.

        num_samples_var: int
            Number of samples to be used in determining signam and mu.
            Only used when approx_mean_var is True.

        Note: init_temp and final_temp are automatically generated if either
              one of the two is not specified.
        """
        if not approx_mean_var:
            # check to see if there are files containing mu and sigma values
            if not os.path.isfile('probe_structure-sigma_mu.npz'):
                self._generate_sigma_mu(num_samples_var)

        print("Generate {} probe structures (struct_per_gen={}, {} present)."
              .format(self.num_to_gen, self.struct_per_gen, self.num_in_gen))
        num_attempt = 0
        while True:
            if atoms is not None:
                self.setting.set_active_template(atoms=atoms,
                                                 generate_template=True)
            else:
                self.setting.set_active_template(size=size,
                                                 unit_cell_id=unit_cell_id,
                                                 generate_template=True)
            # Break out of the loop if reached struct_per_gen
            num_struct = len([row.id for row in
                              self.db.select(gen=self.gen)])
            if num_struct >= self.struct_per_gen:
                break

            struct = self._get_struct_at_conc(conc_type='random')

            print('Generating {} out of {} structures.'
                  .format(num_struct + 1, self.struct_per_gen))
            ps = ProbeStructure(self.setting, struct, self.struct_per_gen,
                                init_temp, final_temp, num_temp,
                                num_steps_per_temp, approx_mean_var)
            probe_struct, cf = ps.generate()
            formula_unit = self._get_formula_unit(probe_struct)
            if self._exists_in_db(probe_struct, formula_unit):
                msg = 'generated structure is already in DB.\n'
                msg += 'generating again... '
                msg += '{} out of {} attempts'.format(num_attempt+1,
                                                      max_attempt)
                print(msg)
                num_attempt += 1
                if num_attempt >= max_attempt:
                    raise MaxAttemptReachedError("Could not generate Emin "
                                                 "structure in {} attempts."
                                                 .format(max_attempt))
            else:
                num_attempt = 0

            kvp = self._get_kvp(probe_struct, cf, formula_unit)
            probe_struct.calc.results['energy'] = None
            self.db.write(probe_struct, kvp)

            if num_attempt >= max_attempt:
                raise MaxAttemptReachedError("Could not generate probe "
                                             "structure in {} attempts."
                                             .format(max_attempt))

    def generate_Emin_structure(self, atoms=None, init_temp=2000,
                                final_temp=1, num_temp=10,
                                num_steps_per_temp=1000,
                                cluster_name_eci=None,
                                random_composition=False):
        """Generate Emin structure.

        Arguments:
        =========
        atoms: Atoms object or a list of Atoms object
            Atoms object with the desired size and composition of the new
            structure. A list of Atoms with different size and/or compositions
            can be passed. Compositions of the supplied Atoms object(s) are
            ignored when random_composition=True.

        init_temp: int or float
            initial temperature (does not represent *physical* temperature)

        final_temp: int or float
            final temperature (does not represent *physical* temperature)

        num_temp: int
            number of temperatures to be used in simulated annealing

        num_steps_per_temp: int
            number of steps in simulated annealing

        cluster_name_eci: dict of list of tuples
            cluster names and ECI values for calculating the energy

        random_composition: bool
            -*False* and atoms = Atoms object: One Emin structure with a
                matching size and composition of the supplied Atoms object is
                generated
            -*False* and atoms = list: The same number of Emin structures that
                matches the length of the list is generated
                Note 1: num_struct_per_gen is ignored and all of the generated
                        structures have the same generation number
                Note 2: each Emin structure will have matching size and
                        composition of the suplied Atoms objects
            -*True* and atoms = Atoms object: Emin structure(s) with a
                matching size of the Atoms object is generated at a random
                composition (within the composition range specified in
                Concentration class)
                Note 1: This will generate Emin structures until the number of
                        structures with the current generation number equals
                        num_struct_per_gen
                Note 2: A check is performed to ensure that none of the newly
                        generated Emin structures have the same composition
            -*True* and atoms = list: The same number of Emin structures that
                matches the length of the list is generated
                Note 1: num_struct_per_gen is ignored and all of the generated
                        structures have the same generation number
                Note 2: each Emin structure will have matching sizes of the
                        supplied Atoms objects but with a random composition
                Note 3: No check is performed to ensure that all new Emin
                        structures have unique composition
        """
        structs = self._set_initial_structures(atoms, random_composition)
        current_count = 0
        num_attempt = 0
        print(structs)
        while current_count < len(structs):
            struct = structs[current_count].copy()
            self.setting.set_active_template(atoms=struct,
                                             generate_template=False)
            print("Generating {} out of {} structures."
                  .format(current_count+1, len(structs)))
            es = EminStructure(self.setting, struct, self.struct_per_gen,
                               init_temp, final_temp, num_temp,
                               num_steps_per_temp, cluster_name_eci)
            emin_struct, cf = es.generate()
            formula_unit = self._get_formula_unit(emin_struct)

            if self._exists_in_db(emin_struct, formula_unit):
                msg = 'generated structure is already in DB.\n'
                msg += 'generating again... '
                msg += '{} out of {} attempts'.format(num_attempt+1,
                                                      max_attempt)
                print(msg)
                num_attempt += 1
                if num_attempt >= max_attempt:
                    raise MaxAttemptReachedError("Could not generate Emin "
                                                 "structure in {} attempts."
                                                 .format(max_attempt))
                continue
            else:
                num_attempt = 0

            print('Structure with E = {:.3f} generated.'.format(es.min_energy))
            kvp = self._get_kvp(emin_struct, cf, formula_unit)
            self.db.write(emin_struct, kvp)

            current_count += 1

    def _set_initial_structures(self, atoms, random_composition=False):
        structs = []
        if isinstance(atoms, Atoms):
            struct = wrap_and_sort_by_position(atoms)
            if random_composition is False:
                self.num_to_gen = 1
                print("Generate 1 Emin structure.")
                structs.append(struct)
            else:
                print("Generate {} Emin structures "
                      "(struct_per_gen={}, {} present)."
                      .format(self.num_to_gen, self.struct_per_gen,
                              self.num_in_gen))
                self.setting.set_active_template(atoms=struct,
                                                 generate_template=True)
                concs = []
                # Get unique concentrations
                num_attempt = 0
                nib = [len(x) for x in self.setting.index_by_basis]
                while len(concs) < self.num_to_gen:
                    x = self.setting.concentration.get_random_concentration(nib=nib)
                    if True in [np.allclose(x, i) for i in concs]:
                        num_attempt += 1
                    else:
                        concs.append(x)
                        num_attempt = 0

                    if num_attempt > 100:
                        raise RuntimeError("Could not find {} unique "
                                           "compositions using the provided "
                                           "Atoms object"
                                           .format(self.num_to_gen))
                num_atoms_in_basis = [len(indices) for indices
                                      in self.setting.index_by_basis]
                for x in concs:
                    num_insert = self.setting.concentration.conc_in_int(
                        num_atoms_in_basis, x)
                    structs.append(self._random_struct_at_conc(num_insert))

        elif all(isinstance(a, Atoms) for a in atoms):
            print("Generate {} Emin structures ".format(len(atoms)))
            if random_composition is False:
                for struct in atoms:
                    structs.append(wrap_and_sort_by_position(struct))
            else:
                concs = []
                nib = [len(x) for x in self.setting.index_by_basis]
                for struct in atoms:
                    self.setting.set_active_template(atoms=struct,
                                                     generate_template=True)
                    x = self.setting.concentration.get_random_concentration(nib=nib)
                    num_atoms_in_basis = [len(indices) for indices
                                          in self.setting.index_by_basis]
                    num_insert = self.setting.concentration.conc_in_int(
                        num_atoms_in_basis, x)
                    structs.append(self._random_struct_at_conc(num_insert))

        else:
            raise ValueError("atoms must be either an Atoms object or a list "
                             "of Atoms objects")
        return structs

    def generate_initial_pool(self):
        """Generate initial pool of random structures."""
        from itertools import product
        print("Generating initial pool consisting of one structure per "
              "concentration where the number of an element is at max/min")
        indx_in_each_basis = []
        start = 0
        for basis in self.setting.concentration.basis_elements:
            indx_in_each_basis.append(list(range(start, start+len(basis))))
            start += len(basis)

        for indx in product(*indx_in_each_basis):
            atoms = self._get_struct_at_conc(conc_type="max", index=indx)
            atoms = wrap_and_sort_by_position(atoms)
            formula_unit = self._get_formula_unit(atoms)

            if not self._exists_in_db(atoms, formula_unit):
                kvp = self.corrfunc.get_cf(atoms)
                kvp = self._get_kvp(atoms, kvp, formula_unit)
                self.db.write(atoms, kvp)

    def _get_struct_at_conc(self, conc_type='random', index=0):
        """Generate a structure at a concentration specified.

        Arguments:
        =========
        conc_type: str
            One of 'min', 'max' and 'random'

        index: int
            index of the flattened basis_element array to specify which element
            to be maximized/minimized
        """
        conc = self.setting.concentration
        if conc_type == 'min':
            x = conc.get_conc_min_component(index)
        elif conc_type == 'max':
            x = conc.get_conc_max_component(index)
        else:
            nib = [len(x) for x in self.setting.index_by_basis]
            x = conc.get_random_concentration(nib=nib)

        num_atoms_in_basis = [len(indices) for indices
                              in self.setting.index_by_basis]
        num_atoms_to_insert = conc.conc_in_int(num_atoms_in_basis, x)
        atoms = self._random_struct_at_conc(num_atoms_to_insert)

        return atoms

    def insert_structure(self, init_struct=None, final_struct=None, name=None,
                         generate_template=False):
        """Insert a user-supplied structure to the database.

        Arguments:
        =========
        init_struct: .xyz, .cif or .traj file
            *Unrelaxed* initial structure.

        final_struct: .traj file (optional)
            Final structure that contains the energy.
            Needs to also supply init_struct in order to use the final_struct.

        name: str (optional)
            Name of the DB entry if a custom name is to be used.
            If *None*, default naming convention will be used.

        generate_template: bool (optional)
            If set to *True*, a template matching the size of the passed
            *init_struct* is created in DB.
        """
        if init_struct is None:
            raise TypeError('init_struct must be provided')

        if name is not None:
            num = sum(1 for _ in self.db.select(['name', '=', name]))
            if num > 0:
                raise ValueError("Name: {} already exists in DB!".format(name))

        if isinstance(init_struct, Atoms):
            init = wrap_and_sort_by_position(init_struct)
        else:
            init = wrap_and_sort_by_position(read(init_struct))

        self.setting.set_active_template(atoms=init_struct,
                                         generate_template=generate_template)

        formula_unit = self._get_formula_unit(init)
        if self._exists_in_db(init, formula_unit):
            raise RuntimeError('Supplied structure already exists in DB')

        cf = self.corrfunc.get_cf(init)
        kvp = self._get_kvp(init, cf, formula_unit)

        if name is not None:
            kvp['name'] = name

        kvp['converged'] = False
        kvp['started'] = False
        kvp['queued'] = False
        kvp['struct_type'] = 'initial'
        uid_init = self.db.write(init, key_value_pairs=kvp)

        if final_struct is not None:
            if isinstance(final_struct, Atoms):
                final = final_struct
            else:
                final = read(final_struct)
            kvp_final = {'struct_type': 'final', 'name': kvp['name']}
            uid = self.db.write(final, kvp_final)
            self.db.update(uid_init, converged=True, started='', queued='',
                           final_struct_id=uid)

    def _exists_in_db(self, atoms, formula_unit=None):
        """Check to see if the passed atoms already exists in DB.

        To reduce the number of assessments for symmetry equivalence,
        check is only performed with the entries with the same concentration
        value.

        Return *True* if there is a symmetry-equivalent structure in DB,
        return *False* otherwise.

        Arguments:
        =========
        atoms: Atoms object

        formula_unit: str
            reduced formula unit of the passed Atoms object
        """
        cond = []
        if formula_unit is not None:
            cond = [("formula_unit", "=", formula_unit)]

        to_prim = True
        try:
            __import__('spglib')
        except Exception:
            msg = "Warning! Setting to_primitive=False because spglib "
            msg += "is missing!"
            print(msg)
            to_prim = False

        symmcheck = SymmetryEquivalenceCheck(angle_tol=1.0, ltol=0.05,
                                             stol=0.05, scale_volume=True,
                                             to_primitive=to_prim)
        atoms_in_db = []
        for row in self.db.select(cond):
            atoms_in_db.append(row.toatoms())
        return symmcheck.compare(atoms.copy(), atoms_in_db)

    def _get_kvp(self, atoms, kvp, formula_unit=None):
        """Get key-value pairs of the passed Atoms object.

        Append terms (started, gen, converged, started, queued, name, conc)
        to key-value pairs and return it.

        Arguments:
        =========
        atoms: Atoms object

        kvp: dict
            key-value pairs (correlation functions of the passed atoms)

        formula_unit: str
            reduced formula unit of the passed Atoms object
        """
        if formula_unit is None:
            raise ValueError("Formula unit not specified!")
        kvp['gen'] = self.gen
        kvp['converged'] = False
        kvp['started'] = False
        kvp['queued'] = False

        count = 0
        for _ in self.db.select(formula_unit=formula_unit):
            count += 1
        kvp['name'] = formula_unit+"_{}".format(count)
        kvp['formula_unit'] = formula_unit
        kvp['struct_type'] = 'initial'
        kvp['unit_cell_id'] = self.setting.unit_cell_id
        kvp['size'] = "x".join((str(d) for d in self.setting.size))
        return kvp

    def _get_formula_unit(self, atoms):
        """Generates a reduced formula unit for the structure."""
        atom_count = []
        all_nums = []
        for group in self.setting.index_by_basis:
            new_count = {}
            for indx in group:
                symbol = atoms[indx].symbol
                if symbol not in new_count.keys():
                    new_count[symbol] = 1
                else:
                    new_count[symbol] += 1
            atom_count.append(new_count)
            all_nums += [v for k, v in new_count.items()]
        gcdp = reduce(lambda x, y: gcd(x, y), all_nums)
        fu = ""
        for i, count in enumerate(atom_count):
            keys = list(count.keys())
            keys.sort()
            for k in keys:
                fu += "{}{}".format(k, int(count[k]/gcdp))
            if i < len(atom_count)-1:
                fu += "_"
        return fu

    def _random_struct_at_conc(self, num_atoms_to_insert):
        """Generate a random structure."""
        rnd_indices = []
        for indices in self.setting.index_by_basis:
            rnd_indices.append(deepcopy(indices))
            shuffle(rnd_indices[-1])

        # Insert the number of atoms
        basis_elem = self.setting.concentration.basis_elements
        assert len(rnd_indices) == len(basis_elem)
        atoms = self.setting.atoms.copy()
        current_conc = 0
        num_atoms_inserted = 0
        for basis in range(len(rnd_indices)):
            current_indx = 0
            for symb in basis_elem[basis]:
                for _ in range(num_atoms_to_insert[current_conc]):
                    atoms[rnd_indices[basis][current_indx]].symbol = symb
                    current_indx += 1
                    num_atoms_inserted += 1
                current_conc += 1
        assert num_atoms_inserted == len(atoms)
        return atoms

    def _determine_gen_number(self):
        """Determine generation number based on the values in DB."""
        try:
            gens = [row.get('gen') for row in self.db.select()]
            gens = [i for i in gens if i is not None]
            gen = max(gens)
            num_in_gen = len([row.id for row in self.db.select(gen=gen)])
            if num_in_gen >= self.struct_per_gen:
                gen += 1
                num_in_gen = 0
        except ValueError:
            gen = 0
            num_in_gen = 0
        return gen, num_in_gen

    def _generate_sigma_mu(self, num_samples_var):
        print('===========================================================\n'
              'Determining sigma and mu value for assessing mean variance.\n'
              'May take a long time depending on the number of samples \n'
              'specified in the *num_samples_var* argument.\n'
              '===========================================================')
        count = 0
        cfm = np.zeros((num_samples_var, len(self.setting.cluster_names)),
                       dtype=float)
        while count < num_samples_var:
            atoms = self._get_struct_at_conc(conc_type='random')
            cfm[count] = self.corrfunc.get_cf(atoms, 'array')
            count += 1
            print('sampling {} ouf of {}'.format(count, num_samples_var))

        sigma = np.cov(cfm.T)
        mu = np.mean(cfm, axis=0)
        np.savez('probe_structure-sigma_mu.npz', sigma=sigma, mu=mu)
