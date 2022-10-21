"""Definition of the Base-class for all Monte Carlo and Simulated Annealing."""
import sys
import numpy as np
from ase.utils import basestring
from ase.units import kB
from ase.montecarlo.swap_atoms import SwapAtoms
from ase.atoms import Atoms
from ase.clease import CEBulk, CECrystal


class MonteCarlo(object):
    """Base-class for all Monte Carlo and Simulated Annealing.

    Arguments
    =========
    atoms: Atoms object to specify the initial structure. A calculator need to
           attached to *atoms* in order to calculate energy.

    setting: None or Setting object for Cluster Expansion.

    temp: temperature in Kelvin for Monte Carlo simulation

    constraint: str
        type of constraints imposed on swapping two atoms.
        - None: any two atoms can be swapped
        - 'nn': any atom selected, swapped only with its nearest neighbor
        - 'basis': any atom selected, swapped only with another atom in the
                   same basis
        - 'nn-basis': any atom selected, swapped only with its nearest neighbor
                      in the same basis

    logfile: file object, str or None.
        - None: logging is disabled
        - string: a file with that name will be opened. If '-', stdout used.
        - file object: use the file object for logging
    """

    def __init__(self, atoms, setting, temp=None, constraint=None,
                 logfile=None):
        if not isinstance(atoms, Atoms):
            raise TypeError('Passed argument should be Atoms object')
        self.atoms = atoms

        if not isinstance(setting, (CEBulk, CECrystal, None)):
            raise TypeError("setting must be CEBulk or CECrystal "
                            "object for Cluster Expansion. Set as *None* "
                            "otherwise.")
        self.setting = setting
        self.energy = None

        if isinstance(temp, (int, float)):
            self.kT = kB * temp
        elif temp is None:
            self.kT = None
        else:
            raise TypeError('temp needs to be int of float type')

        # constraints
        allowed_constraints = [None, 'nn', 'basis', 'nn-basis']
        if constraint not in allowed_constraints:
            raise TypeError('constraint needs to be one of: '
                            '{}'.format(allowed_constraints))
        self.constraint = constraint

        # logfile
        if isinstance(logfile, basestring):
            if logfile == '-':
                logfile = sys.stdout
            else:
                logfile = open(logfile, 'a')
        self.logfile = logfile
        self.nsteps = 0

    def run(self, num_steps=100, average=False):
        """Run Monte Carlo / Simulated Annealing.

        Arguments
        =========
        num_steps: Number of steps in Monte Carlo simulation.
        average: whether or not to return the average values.
            - True: returns average energy and energy**2 over entire simulation
            - False: returns the sum of energy and energy**2 over the entire
                     simulation
        """
        # starting energy
        self.energy = self.atoms.get_potential_energy()

    def _swap(self):
        """Swap two atoms.

        Evaluate whether or not to accept the new structure. If not accepted,
        revert back to the orignal structure. Regardless of the acceptance,
        return the energy of the swapped structure
        """
        swap = SwapAtoms(self.setting)
        if self.constraint is None:
            swapped_indices = swap.swap_any_two_atoms(self.atoms)
        elif self.constraint == 'nn':
            swapped_indices = swap.swap_nn_atoms(self.atoms)
        elif self.constraint == 'basis':
            swapped_indices = swap.swap_any_two_in_same_basis(self.atoms)
        elif self.constraint == 'nn-basis':
            swapped_indices = swap.swap_nn_in_same_basis(self.atoms)
        else:
            raise NotImplementedError('This feature is not implemented')

        energy = self.atoms.get_potential_energy()
        accept = np.exp((self.energy - energy) / self.kT) > np.random.uniform()

        if accept:
            self.energy = energy
        else:
            # Swap atoms back to the original
            swap.swap_by_indices(self.atoms, swapped_indices[0],
                                 swapped_indices[1])
            # CE calculator needs to call a *restore* method
            if self.atoms.calc.__class__.__name__ == 'Clease':
                self.atoms.calc.restore()

        self.nsteps += 1

        return accept, energy

    def log(self, accept=None, new_energy=None, logtype='MC'):
        """Writes log to to either a file or stdout.

        If logtype = 'MC', write accept, energy and energy^2 of candidate and
        selected strucutres.

        Otherwise (i.e., Simulated Annealing), write kT, accept and energy
        """
        if self.logfile is None:
            return True

        if self.nsteps == 0:
            if logtype == 'MC':
                self.logfile.write('\t\t\t\tselected structure \t\t\t'
                                   'candidate structure\n')
                self.logfile.write('steps \taccept \tEnergy \t\t\tEnergy^2'
                                   '\t\tEnergy \t\t\tEnergy^2\n')
                self.logfile.write('{}\t{}\t{}\t{}'.format(self.nsteps,
                                                           "-----",
                                                           self.energy,
                                                           self.energy**2))
            else:
                self.logfile.write('steps \tkT \t\t\taccept \tEnergy\n')
                self.logfile.write('{}\t{}\t\t\t{}\t{}'.format(self.nsteps,
                                                               '-----',
                                                               '-----',
                                                               self.energy))

        else:
            if logtype == 'MC':
                self.logfile.write('{}\t{}\t{}'.format(self.nsteps, accept,
                                                       self.energy))
                self.logfile.write('\t{}\t{}\t{}'.format(self.energy**2,
                                                         new_energy,
                                                         new_energy**2))
            else:
                self.logfile.write('{}\t{}\t{}\t{}'.format(self.nsteps,
                                                           self.kT, accept,
                                                           self.energy))
        self.logfile.write('\n')
        self.logfile.flush()
