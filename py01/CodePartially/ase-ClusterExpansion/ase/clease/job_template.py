import os
from ase.db import connect
from ase.io import write
from subprocess import check_output
from ase.clease import jobscript_template as js


class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


class Submit(object):
    """
    Class that handles submission of DFT computations to the batch system.
    """
    def __init__(self, db_name, compute_param1, compute_param2, num_cores,
                 num_nodes, email):
        self.db_name = db_name
        self.compute_param1 = compute_param1
        self.compute_param2 = compute_param2
        self.num_cores = int(num_cores)
        self.num_nodes = int(num_nodes)
        self.email = email
        self.db = connect(self.db_name)

    @property
    def submit_new(self):
        """
        Find a row in the database that has the condition 'queued=False' and
        'started=False' and submit it to a SLURM batch system.
        """
        db_name = self.db_name
        ids = [row.id for row in
               self.db.select(struct_type='initial',
                              queued=False, started=False)]
        id = ids[0]
        name = self.db.get(id).name
        atoms = self.db.get_atoms(id)
        self.compute_param1['system'] = name
        self.compute_param2['system'] = name
        suffix_length = len(name.split('_')[-1]) + 1
        prefix = name[:(-1*suffix_length)]
        self.remove_X(atoms)

        if not os.path.exists('./{}/{}'.format(prefix, name)):
            os.makedirs('./{}/{}'.format(prefix, name))
            os.system('chmod -R 755 ./{}/{}'.format(prefix, name))

        with cd('./{}/{}'.format(prefix, name)):
            # write traj file
            write('input.traj', atoms)
            vasp_file = open('vasp.py', 'w')
            vasp_file.write(js.vasp_new(self.compute_param1,
                                        self.compute_param2,
                                        name, db_name))
            vasp_file.close()

            runfile = open('run.sh', 'w')
            if self.num_cores == 8:
                runfile.write(js.slurm_script_8(name, self.num_nodes,
                                                self.email))
            elif self.num_cores == 16:
                runfile.write(js.slurm_script_16(name, self.num_nodes,
                                                 self.email))
            else:
                runfile.write(js.slurm_script_24(name, self.num_nodes,
                                                 self.email))
            runfile.close()

            os.system('chmod u+xwr vasp.py run.sh')
            output_string = check_output(['sbatch run.sh'], shell=True,
                                         universal_newlines=True)

        # wait and get the response to confirm
        if "Submitted" not in output_string:
            raise ValueError("Job name {} may not be submitted".format(name))

        print('submitted {}'.format(name))
        self.db.update(id, queued=True)
        return True

    @property
    def submit_restart(self):
        db_name = self.db_name
        ids = [row.id for row in
               self.db.select(struct_type='initial',
                              converged=False, started=True)]
        all_names = []
        for id in ids:
            all_names.append(self.db.get(id).name)

        queue_list = self.jobs_in_queue
        names = [i for i in all_names if i not in queue_list]
        indices = [all_names.index(name) for name in names]
        id = ids[indices[0]]
        name = names[0]
        self.compute_param1['system'] = name
        self.compute_param2['system'] = name

        print('There are {} jobs to be restarted.'.format(len(names)))
        print('Names of jobs are:')
        print(names)

        suffix_length = len(name.split('_')[-1]) + 1
        prefix = name[:(-1*suffix_length)]
        with cd('./{}/{}'.format(prefix, name)):
            # Check the reason for incompletion
            if os.path.isfile('{}.log'.format(name)):
                logfile = open('{}.log'.format(name), 'r')
                log_msg = logfile.read()
                logfile.close()
                # Case where the last job halted due to the time limit
                if "TIME LIMIT" in log_msg:
                    print('Job previously halted due to the walltime limit.')
                    print('Resuming with regular conditions.')
                # Case where there was some error in the previous run
                else:
                    raise ValueError('This job was halted due to an error.\n'
                                     'Delete the .log file to continue after '
                                     'the issues are fixed.')

            vasp_file = open('vasp.py', 'w')
            vasp_file.write(js.vasp_restart(self.compute_param1,
                                            self.compute_param2,
                                            name, db_name))
            vasp_file.close()

            runfile = open('run.sh', 'w')
            if self.num_cores == 8:
                runfile.write(js.slurm_script_8(name, self.num_nodes,
                                                self.email))
            elif self.num_cores == 16:
                runfile.write(js.slurm_script_16(name, self.num_nodes,
                                                 self.email))
            else:
                runfile.write(js.slurm_script_24(name, self.num_nodes,
                                                 self.email))
            runfile.close()

            os.system('chmod u+xwr vasp.py run.sh')
            output_string = check_output(['sbatch run.sh'], shell=True,
                                         universal_newlines=True)

        # wait and get the response to confirm
        if "Submitted" not in output_string:
            raise ValueError("Job name {} may not be submitted".format(name))

        print('submitted {}'.format(name))
        self.db.update(id, queued=True, converged=False, started=False)
        return True

    def remove_X(self, atoms):
        """
        Vacancies are specified with the ghost atom 'X', which must be removed
        before passing the atoms object to a calculator.
        """
        del atoms[[atom.index for atom in atoms if atom.symbol == 'X']]

    @property
    def jobs_in_queue(self):
        """
        Returns a list of job names that are in the SLURM batch system.
        """
        job_names = []
        jobs_string = check_output(['qstat -f | grep -C 1 $USER'],
                                   shell=True,
                                   universal_newlines=True).splitlines()
        for line in jobs_string:
            if 'Job_Name' not in line:
                continue
            line = line.split()
            job_names.append(line[-1])
        return job_names
