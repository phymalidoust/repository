#!/bin/sh
#SBATCH --partition=WORKQ
#SBATCH --time=4-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name="hcpAMZ"
#SBATCH --output="output8.out"

module purge
module load GCC/6.3.0-2.27  
module load OpenMPI/2.0.2
module load icc/2017.1.132-GCC-6.3.0-2.27
module load impi/2017.1.132
module load ifort/2017.1.132-GCC-6.3.0-2.27
module load impi/2017.1.132
module load Python/3.6.1

export PATH=${PATH}:"/home/mohamal/MPI/mpich3-install/bin"
export PATH=${PATH}:"/home/mohamal/python/bin"
export PYTHONPATH=${PYTHONPATH}:"/home/mohamal/python/Python-3.6.0"
export PATH=${PATH}:"/home/mohamal/.local/bin"
export PYTHONPATH=${PYTHONPATH}:"/home/mohamal/ase"
export PATH=${PATH}:"/home/mohamal/ase/bin"
export PATH=${PATH}:"/home/mohamal/gcc-8.2.0/bin"
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:"/home/mohamal/gcc-8.2.0/lib"
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:"/home/mohamal/gcc-8.2.0/lib64"



cd /home/mohamal/RUNs/runs_ce_parall/hcp/AlMgZn

mpirun -np $SLURM_NTASKS python3 /home/mohamal/RUNs/runs_ce_parall/hcp/AlMgZn/parl_3elemnts.py 8 24










