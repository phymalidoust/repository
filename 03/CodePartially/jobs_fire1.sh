#!/bin/bash
#SBATCH --job-name=AMZ
#SBATCH --time=0-24:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --account=nn9497k
#SBATCH --mem=32000MB
module purge
module load FFTW/3.3.6-intel-2016b
module load intel/2016b
module load Python/2.7.12-intel-2016b
export PATH=${PATH}:"/home/davidkl/.local/bin"
export PYTHONPATH=${PYTHONPATH}:"/home/davidkl/.local/lib/python2.7/site-packages"
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:"/home/davidkl/.local/lib"
export GPAW_SETUP_PATH="/home/davidkl/gpaw-setups-0.9.20000"
WORKDIR=/global/work/mohamal
cd ${WORKDIR}
mpirun -np $SLURM_NTASKS gpaw-python /home/mohamal/gpaw_test/hcp/AlMgZn/dft_fire1.py 145 2.8 fire /home/mohamal/gpaw_test/hcp/AlMgZn/almgzn_hcp_2x2x2.db
