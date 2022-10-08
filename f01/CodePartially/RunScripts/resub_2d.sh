#!/bin/bash
echo
echo The number of nodes and cores you supplied will not work
echo for  the given 2D problem conditions
echo Trying to find a configuration that works
echo

MAX_CORES_NODE=$(( BC_STANDARD_NODE_CORES / 2 ))

FINISHED=false
for ((TRY_NODES = $NODES; TRY_NODES <= $((NODES + 20)); TRY_NODES++ )); do
    for (( SET_NODE_CORES=4; SET_NODE_CORES<=$MAX_CORES_NODE; SET_NODE_CORES++ )); do

	# Most of the HPC systems only allow MPI procs to evenly divide
	# the total cores on a node so enforce that condition
	if [ $BC_HOST == 'onyx' ]; then
	    CORES_MOD=$(( BC_STANDARD_NODE_CORES % SET_NODE_CORES ))
	    if [ $CORES_MOD -ne 0 ]; then
		continue
	    fi	    
	fi

	TOTAL_CORES=$(( TRY_NODES * SET_NODE_CORES ))

        ans=$(python3 launch_params.py $TRY_NODES $TOTAL_CORES | grep omp_threads)
	THREADS=$(echo $ans | grep omp_threads | cut -d : -f 2)

	if [ $THREADS != 100 ]; then
	    if [ -e resub_2d.pbs ]; then
		rm resub_2d.pbs
	    fi

	    # Need a way to determine a reasonable new wall time for the submitted job
	    WALLTIME=$(qstat -f $PBS_JOBID | sed -rn 's/.*Resource_List.walltime = (.*)/\1/p')	
	    if [ $PBS_QUEUE == "debug" ]; then
		WALLTIME="1:00:00"
	    fi

	    cat >> resub_2d.pbs <<EOF 
#!/bin/bash
#PBS -l select=${TRY_NODES}:ncpus=${BC_STANDARD_NODE_CORES}:mpiprocs=${SET_NODE_CORES}
#PBS -l walltime=${WALLTIME}
#PBS -q ${PBS_QUEUE}
#PBS -j oe
#PBS -N spin_flip
#PBS -A ${PBS_ACCOUNT}
#**************************************************************************
# USER-SPECIFIED INPUTS
#**************************************************************************

# Set desired workdir destination
# This script will create a job-specific sub-directory in this location
# which will store all output data
workpath=$WORKDIR/Spintronics

# Would you like to restart from g.txt from a previous run?
export spin_restart=${spin_restart}
# Magnetic field for each ferromagnet in order
export spin_H0="${spin_H0}"
# Minimum number of x-momenta to solve (= # of matrices to solve)
export spin_min_nkx=0
# Minimum number kx factor
export spin_min_kx_factor=${spin_min_kx_factor}
# Max iterations to attempt if relative error doesn't decrease enough
export spin_itmax=${spin_itmax}
#Convergence tolerance
export spin_minerror=${spin_minerror}
# Time
export spin_time=${spin_time}
# Non-default relative temperature
export spin_rel_temp=${spin_rel_temp}
# Dimensions of the domain in nondimensional length
export spin_kfl1=${spin_kfl1}
export spin_kfl2=${spin_kfl2}
export spin_kfl3=${spin_kfl3}
export spin_kfds=${spin_kfds}
export spin_kfw=${spin_kfw}
# Phase of the Josephson Junction with in the superconductors
export spin_deg_vphi2=${spin_deg_vphi2}
export spin_deg_vphi1=${spin_deg_vphi1}
export spin_vphi1=${spin_vphi1}
export spin_vphi2=${spin_vphi2}
# Magnetic field exchange angles
export spin_deg_theta="${spin_deg_theta}"
export spin_deg_phi="${spin_deg_phi}"
export spin_phi="${spin_phi}"
export spin_theta="${spin_theta}"

# Input the range of energies to output to W_.txt files
export spin_Wtxt_range=${spin_Wtxt_range}

# Impurity scattering strength (H_B)
export spin_scatt_str=${spin_scatt_str}
# Impurity concentrations for each section [S1,F1,F2,F3,S2]
export spin_imp_conc=${spin_imp_conc}
# Option for seed for random number generator used to
# determine impurity distribution
# 0: default const. seed (always repeatable)
# 1: inputted const. seed (user-controlled repeatable)
# 2: random seed (never repeatable)
export spin_rand_opt=${spin_rand_opt}
# Seed for random number generator with option 1
export spin_rand_seed=${spin_rand_seed}

# Range of energies for DOS calculation
export spin_DOS_range="${spin_DOS_range}"
# Interval between energies for DOS
export spin_DOS_inter=${spin_DOS_inter}
# Three energy states for DOS output
export spin_DOS_energy="${spin_DOS_energy}"
# Three (y,z) locations for DOS output ordered [y1,y2,y3] & [z1,z2,z3]
y1=\$(python -c "print(\$spin_kfds/2.0)")
y2=\$(python -c "print(\$spin_kfds+\$spin_kfl1/2.0)")
y3=\$(python -c "print(\$spin_kfds+\$spin_kfl1+\$spin_kfl2/2.0)")
export spin_DOS_loc_y="[\$y1,\$y2,\$y3]"

z1=\$(python -c "print(\$spin_kfw/2.0)")
z2=\$(python -c "print(\$spin_kfw/2.0)")
z3=\$(python -c "print(\$spin_kfw/2.0)")
export spin_DOS_loc_z="[\$z1,\$z2,\$z3]"


#**************************************************************************
# SYSTEM AND JOB SETUP
#**************************************************************************
cd \$PBS_O_WORKDIR
SRC=\$PBS_O_WORKDIR
exe=spin_flip_2d

# Call external script that loads system dependent setup 
source sys_setup.sh

echo Nodes = \$NODES

export NSLOTS=\$(wc -l \$PBS_NODEFILE | cut -f1 -d " ")
echo Cores = \$NSLOTS

ans=\$(python3 \${SRC}/launch_params.py \$NODES \$NSLOTS | grep omp_threads)
THREADS=\$(echo \$ans | grep omp_threads | cut -d : -f 2)

let "RANKS_PER_NODE = \$NSLOTS / \$NODES"
let "RANKS = \$RANKS_PER_NODE * \$NODES"
export OMP_NUM_THREADS=\$THREADS

echo Threads = \$THREADS
echo Ranks per Node = \$RANKS_PER_NODE
echo Ranks = \$RANKS
echo

# create a folder for this job
job_num=\$(echo \$PBS_JOBID | cut -d '.' -f 1)
job_dir=\$workpath/\$job_num
mkdir -p \$job_dir

cp \$PBS_O_WORKDIR/\${exe} \$job_dir/
cp \$PBS_O_WORKDIR/input.txt \$job_dir/
cp \$PBS_O_WORKDIR/generate_report.py \$job_dir/

# link to the most recent job
rm -rf \$WORKDIR/Spintronics/latest
ln -s \$job_dir \$WORKDIR/Spintronics/latest

cd \$job_dir

ulimit -c unlimited

#**************************************************************************
# RUN SPIN_FLIP_2D
#**************************************************************************

if check_module PrgEnv-gnu; then
   aprun -d \$THREADS -n \$RANKS -N \$RANKS_PER_NODE  \${job_dir}/\${exe}
else
   mpiexec_mpt -n \$RANKS \${job_dir}/\${exe}
fi

python3 generate_report.py


EOF
           echo
	   echo Creating new PBS script with $TRY_NODES nodes and $SET_NODE_CORES cores per node
	   FINISHED=true
	   break
	fi # THREADS != 100
    done   # Cores
    if $FINISHED ; then
	break
    fi
done       # Nodes

if [ $THREADS == 100 ]; then
    echo
    echo Cannot find a number of cores that works for this 2D problem
    echo
fi
