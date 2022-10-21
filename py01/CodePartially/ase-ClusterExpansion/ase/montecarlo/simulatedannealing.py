"""Definition of SimulatedAnnealing Class."""
import math
import numpy as np
from ase.units import kB
from ase.montecarlo.montecarlo import MonteCarlo


class SimulatedAnnealing(MonteCarlo):
    """Class for performing Simulated Annealing.

    Arguments
    =========
    atoms: Atoms object to specify the initial structure. A calculator need to
           attached to *atoms* in order to calculate energy.

    init_temp: initial temperature in Kelvin.

    final_temp: final temperature in Kelvin.

    num_temp: number of temperatures to be used, including init_temp and
              final_temp

    interpolation: interpolation scheme to be used for generating intermediate
                   temptratures
        - 'linear': use linear-scale interpolation
        - 'log': use log-scale interpolation

    temp_list: list or numpy array containing temperatures to be used. If
               temp_list is given, init_temp, final_temp, num_temp and
               interpolation are ignored.

    constraint: types of constraints imposed on swapping two atoms.
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

    def __init__(self, atoms, setting=None, init_temp=None, final_temp=None,
                 num_temp=None, interpolation='linear', temp_list=None,
                 constraint=None, logfile=None):
        MonteCarlo.__init__(self, atoms=atoms, setting=setting, temp=None,
                            constraint=constraint, logfile=logfile)

        if isinstance(temp_list, (list, np.ndarray)):
            # case for using temp_list
            # sanity check
            if len(temp_list) < 1:
                raise ValueError('length of temp_list must be > 0')
            elif not all(isinstance(i, (int, float)) for i in temp_list) \
                    and isinstance(temp_list, list):
                raise TypeError('all items in list need to be int or float')
            elif not all(i > 0 for i in temp_list):
                raise ValueError('all temperature values need to be positive')
            elif not all(earlier >= later for earlier, later in
                         zip(temp_list, temp_list[1:])):
                raise ValueError('temperatures must be in descending order')
            # set temperatures
            self.kTs = kB * np.array(temp_list)
            self.num_temp = len(self.kTs)

        elif temp_list is not None:
            raise TypeError('temp_list should be a list, numpy array, or None '
                            'type')

        else:
            # case for using init_temp, final_temp, num_temp
            # sanity check
            if not isinstance(init_temp, (int, float)) or \
               not isinstance(final_temp, (int, float)):
                raise TypeError('temperatures must be int or float type')
            if init_temp <= final_temp:
                raise ValueError('init_temp must be higher than final_temp')
            if not isinstance(num_temp, int):
                raise TypeError('num_temp must be int type')
            if not (init_temp > 0 and final_temp > 0 and num_temp > 0):
                raise ValueError('init_temp, final_temp and num_temp must be '
                                 'positive')

            # set temperatures
            if interpolation == 'linear':
                self.kTs = kB * np.linspace(init_temp, final_temp, num_temp)
            elif interpolation == 'log':
                self.kTs = kB * np.logspace(math.log10(init_temp),
                                            math.log10(final_temp),
                                            num_temp)
            else:
                raise ValueError("interpolation must be 'linear' or 'log'")
            self.num_temp = num_temp

    def run(self, num_steps=100):
        """Run Simulated Annealing.

        Returns the minimum energy structure (Atoms object) and its energy.

        Arguments
        =========
        num_steps: Number of steps in Simulated Annealing.
        """
        # starting energy
        MonteCarlo.run(self)
        self.log(logtype='SA')

        if num_steps < self.num_temp:
            raise ValueError('number of Simulated Annealing steps must be '
                             'larger than the number of temperatures')
        steps_per_temp = int(num_steps / self.num_temp)

        temp_indx = 0
        self.kT = self.kTs[temp_indx]

        while self.nsteps < num_steps:
            accept, energy = self._swap()
            self.log(accept, energy, 'SA')
            if self.nsteps % steps_per_temp == 0 and \
               temp_indx < self.num_temp - 1:
                temp_indx += 1
                self.kT = self.kTs[temp_indx]
            # keep track of latest accepted structure
            if accept:
                atoms = self.atoms.copy()
                final_energy = energy

        return atoms, final_energy
