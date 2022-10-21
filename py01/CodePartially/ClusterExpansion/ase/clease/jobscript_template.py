from textwrap import dedent


def vasp_new(param1, param2, name, db_name):
    script_text = """\
    import os
    from ase.calculators.vasp import Vasp2
    from ase.io import read, write
    from ase.db import connect
    from ase.io.trajectory import TrajectoryWriter
    from ase.clease.tools import update_db

    # update database
    name='"""+str(name)+"""'
    db_name='../../"""+str(db_name)+"""'
    db = connect(db_name)
    uid_initial = db.get(name=name, struct_type='initial').id
    db.update(uid_initial, started=True, queued=False)

    # run calculation
    atoms=read('input.traj')
    compute_param="""+str(param1)+"""
    calc = Vasp2(**compute_param)
    atoms.set_calculator(calc)
    energy = atoms.get_potential_energy()
    # store in traj file
    outcar = read('OUTCAR', index=':')
    writer = TrajectoryWriter('output.traj', mode='w')
    for atoms in outcar:
        writer.write(atoms=atoms)
    os.system('rm WAVECAR')


    # run with the final structure to ensure convergence
    atoms = read('output.traj', -1)
    compute_param="""+str(param2)+"""
    calc = Vasp2(**compute_param)
    atoms.set_calculator(calc)
    energy = atoms.get_potential_energy()
    # append the image to traj file
    outcar = read('OUTCAR', index=':')
    writer = TrajectoryWriter('output.traj', mode='a')
    for atoms in outcar:
        writer.write(atoms=atoms)

    # It is converged at this point -> update database
    update_db(uid_initial=uid_initial,
              final_struct=read('output.traj', -1),
              db_name=db_name)
    """
    return dedent(script_text)


def vasp_restart(param1, param2, name, db_name):
    script_text = """\
    import os
    from ase.calculators.vasp import Vasp2
    from ase.io import read, write
    from ase.db import connect
    from ase.io.trajectory import TrajectoryWriter
    from ase.clease.tools import update_db

    # update database
    name='"""+str(name)+"""'
    db_name='../../"""+str(db_name)+"""'
    db = connect(db_name)
    uid_initial = db.get(name=name, struct_type='initial').id
    db.update(uid_initial, started=True, queued=False)

    # save the OUTCAR from the previous run to traj
    outcar = read('OUTCAR', index=':')
    writer = TrajectoryWriter('output.traj', mode='a')
    for atoms in outcar:
        writer.write(atoms=atoms)

    # run calculation
    atoms = read('output.traj', -1)
    compute_param="""+str(param1)+"""
    calc = Vasp2(**compute_param)
    atoms.set_calculator(calc)
    energy = atoms.get_potential_energy()
    # store in traj file
    outcar = read('OUTCAR', index=':')
    for atoms in outcar:
        writer.write(atoms=atoms)

    # run with the final structure to ensure convergence
    atoms = read('output.traj', -1)
    compute_param="""+str(param2)+"""
    calc = Vasp2(**compute_param)
    atoms.set_calculator(calc)
    energy = atoms.get_potential_energy()
    # append the image to traj file
    outcar = read('OUTCAR', index=':')
    writer = TrajectoryWriter('output.traj', mode='a')
    for atoms in outcar:
        writer.write(atoms=atoms)

    # It is converged at this point -> update database
    update_db(uid_initial=uid_initial,
              final_struct=read('output.traj', -1),
              db_name=db_name)
    """
    return dedent(script_text)


def slurm_script_8(job_name, num_nodes, email):
    script_text = """\
    #!/bin/bash
    #SBATCH --mail-user="""+str(email)+"""
    #SBATCH --mail-type=ALL
    #SBATCH -N """+str(num_nodes)+"""
    #SBATCH -n """+str(int(num_nodes*8))+"""
    #SBATCH --time=150:00:00
    #SBATCH --output="""+str(job_name)+""".log
    #SBATCH --job-name="""+str(job_name)+"""
    #SBATCH --partition=xeon8

    module load VASP/5.4.1-iomkl-2017a
    python vasp.py"""
    return dedent(script_text)


def slurm_script_16(job_name, num_nodes, email):
    script_text = """\
    #!/bin/bash
    #SBATCH --mail-user="""+str(email)+"""
    #SBATCH --mail-type=ALL
    #SBATCH -N """+str(num_nodes)+"""
    #SBATCH -n """+str(int(num_nodes*16))+"""
    #SBATCH --time=150:00:00
    #SBATCH --output="""+str(job_name)+""".log
    #SBATCH --job-name="""+str(job_name)+"""
    #SBATCH --partition=xeon16

    module load VASP/5.4.1-iomkl-2017a
    python vasp.py"""
    return dedent(script_text)


def slurm_script_24(job_name, num_nodes, email):
    script_text = """\
    #!/bin/bash
    #SBATCH --mail-user="""+str(email)+"""
    #SBATCH --mail-type=ALL
    #SBATCH -N """+str(num_nodes)+"""
    #SBATCH -n """+str(int(num_nodes*24))+"""
    #SBATCH --time=50:00:00
    #SBATCH --output="""+str(job_name)+""".log
    #SBATCH --job-name="""+str(job_name)+"""
    #SBATCH --partition=xeon24

    module load VASP/5.4.1-iomkl-2017a
    python vasp.py"""
    return dedent(script_text)
