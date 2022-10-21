"""Definition of Metropolis Class for Metropolis Monte Carlo Simulations."""
from ase.montecarlo.montecarlo import MonteCarlo


class Metropolis(MonteCarlo):
    """Class for performing Metropolis Monte Carlo sampling.

    Arguments
    =========
    atoms: Atoms object to specify the initial structure. A calculator need to
           attached to *atoms* in order to calculate energy.

    setting: None or Setting object for Cluster Expansion.

    temp: temperature in Kelvin for Monte Carlo simulation

    constraint: types of constraints imposed on swapping two atoms.
        - None: any two atoms can be swapped
        - 'nn': any atom selected, swapped only with its nearest neighbor
        - 'basis': any atom selected, swapped only with another atom in the
                   same basis
        - 'nn-basis': any atom selected, swapped only with its nearest neighbor
                      in the same basis

    logfile: file object or str
        If *logfile* is a string, a file with that name will be opened.
        Use '-' for stdout.
    """

    def __init__(self, atoms, setting=None, temp=293, constraint=None,
                 logfile=None):
        if temp is None:
            raise TypeError('temp needs to be int or float type')
        MonteCarlo.__init__(self, atoms=atoms, setting=setting, temp=temp,
                            constraint=constraint, logfile=logfile)

    def run(self, num_steps=10, average=False):
        """Run Monte Carlo simulation.

        Perform Metropolis Monte Carlo simulation using the number of steps
        specified by a user. Returns energy and energy**2.

        Arguments
        =========
        num_steps: Number of steps in Monte Carlo simulation.
        average: whether or not to return the average values.
            - True: returns average energy and energy**2 over entire simulation
            - False: returns the sum of energy and energy**2 over the entire
                     simulation
        """
        # starting energy
        MonteCarlo.run(self)
        energy_sum = self.energy
        energy_sq_sum = self.energy**2
        self.log()

        while self.nsteps < num_steps:
            accept, energy = self._swap()
            energy_sum += self.energy
            energy_sq_sum += self.energy**2
            self.log(accept, energy)

        if average:
            energy_sum /= self.nsteps
            energy_sq_sum /= self.nsteps

        return energy_sum, energy_sq_sum
