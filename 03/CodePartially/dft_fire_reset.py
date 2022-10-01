import gpaw as gp
from ase.db import connect
#from ase.optimize.precon import PreconLBFGS
from ase.optimize.precon import PreconFIRE
from ase.io.trajectory import Trajectory
from ase.calculators.singlepoint import SinglePointCalculator
import sys

folder = "/global/work/mohamal"
class RestartSave(object):
    def __init__(self, name, calc):
        self.calc = calc
        self.fname = folder + '/' + name + '.gpw'

    def __call__(self):
        self.calc.write(self.fname, mode='all')

def main(argv):
    runID = int(argv[0])  # ID in the SQL database
    kpt_density = float(argv[1])

    print("Running job: {}".format(runID))
    db_name = "/home/mohamal/gpaw_test/hcp/AlMgZn/almgzn_hcp_2x2x2.db"

    # Create a database connection
    db = connect(db_name)

    # Retrieve the unique name of this particular entry
    name = db.get(id=runID).key_value_pairs["name"]

    # Update the databse
    db.update(runID, started=True, converged=False)
    db.update(runID, kpt_density=kpt_density)

    # Get the atoms object from the database
    #atoms = db.get_atoms(id=runID)

    #nbands = "120%"  # Number of electronic bands

    #kpts = {"density": kpt_density, "even": True}
    #calc = gp.GPAW(mode=gp.PW(600), xc="PBE", kpts=kpts, nbands=nbands)

    kpts = {"density": kpt_density, "even": True}
    fname = folder + '/' + name + '.gpw'

    atoms, calc = gp.restart(fname, kpts=kpts)
    #calc = gp.GPAW(fname)

    restart_save = RestartSave(name, calc)
    
    atoms.set_calculator(calc)

    logfile = "log_fir{}.log".format(runID)
    traj = "traj_fir{}.traj".format(runID)

    # Initialize a trajactory file
    trajObj = Trajectory(traj, 'w', atoms)

    # Initialize the relaxer
    relaxer = PreconFIRE(atoms, logfile=logfile, use_armijo=True,
                          variable_cell=True)
    relaxer.attach(restart_save)
    # Attach the trajectory file to the relaxer
    relaxer.attach(trajObj)

    # Run until force and stress criteria have been met
    
    fmax = 0.03  # Maximum force in eV/Å
    smax = 0.01  # Maximum stress in eV/Å^3

##    fmax = 0.05  # Maximum force in eV/Å
##    smax = 0.05  # Maximum stress in eV/Å^3

    
    relaxer.run(fmax=fmax, smax=smax)

    # Get and print the total energy
    energy = atoms.get_potential_energy()
    print("Energy: {}".format(energy))

    # What follows is very crucial that it is done exactly like this

    # Retrieve the original (unrelaxed object) from the database
    orig_atoms = db.get_atoms(id=runID)

    # Attacha singlet point calculator with the energy of the relaxed structure
    scalc = SinglePointCalculator(orig_atoms, energy=energy)
    orig_atoms.set_calculator(scalc)

    # Get a all the key_value_pairs
    kvp = db.get(id=runID).key_value_pairs

    # Delete the original entry
    del db[runID]

    # Write the new objet to the database
    # Unrelaxed system, with the energy of the relaxed one
    newID = db.write(orig_atoms, key_value_pairs=kvp)

    # Update the converged flag
    db.update(newID, converged=True)

    # Store also the relaxed object (NOTE: they are linked via their name)
    db.write(atoms, name=name, state="relaxed")


if __name__ == "__main__":
    main(sys.argv[1:])
