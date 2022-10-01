#!/bin/bash

# function for checking if a module is available
check_module(){
    result=$(module show $1 2>&1)
    if [[ $result == *"Unable to locate a modulefile for '$1'"* ]]; then
	return 1    # false
    else
	return 0    # true
    fi
}

#if check_module PrgEnv-cray; then
if [ $BC_HOST == 'onyx' ]; then
    # we're on a cray system
    module unload PrgEnv-cray PrgEnv-gnu PrgEnv-pgi PrgEnv-intel
    module load PrgEnv-gnu
    
    module load cray-libsci
    module load cray-hdf5
    module load cseinit
    module load cse/anaconda3

    #ensure latest version of gcc compiler loaded
    module unload gcc/4.9.3
    module load gcc/9.3.0

    NODES=$(aprun -q -B hostname | sort | uniq | wc -l)

elif [ $BC_HOST == 'narwhal' ]; then

    module unload PrgEnv-cray PrgEnv-intel PrgEnv-gnu
    module load PrgEnv-gnu

    module unload gcc/10.3.0
    module load gcc/9.3.0
    module load cray-hdf5
    module load cray-libsci
    module load python/3
    
    NODES=$(cat $PBS_NODEFILE | cut -d " " -f 1 | sort | uniq | wc -l | cut -d " " -f 1)
    
else
    # We are on an SGI system
    module load cseinit
    module load cse/anaconda3
    module unload compiler/intel gcc
    module unload cse/openmpi mpt
    module load gcc
    module load mpt
    module load costinit
    module load hdf5/gnu
    if [ $BC_HOST == "mustang" ]; then
	MKLROOT=/p/app/intel/parallel_studio_xe_2018_update3/compilers_and_libraries_2018.3.222/linux/mkl/lib/intel64
	INTELROOT=/p/app/intel/parallel_studio_xe_2018_update3/compilers_and_libraries_2018.3.222/linux/compiler/lib/intel64
    elif [ $BC_HOST == "gaffney" ]; then
	MKLROOT=/p/app/intel/parallel_studio_xe_2019_update4/compilers_and_libraries/linux/mkl/lib/intel64
	INTELROOT=/p/app/intel/parallel_studio_xe_2019_update4/compilers_and_libraries/linux/compiler/lib/intel64
    elif [ $BC_HOST == "centennial" ]; then
	MKLROOT=/p/app/intel/l_psxe_2017.1.043/compilers_and_libraries_2017.1.132/linux/mkl/lib/intel64
	INTELROOT=/p/app/intel/l_psxe_2017.1.043/compilers_and_libraries_2017.1.132/linux/compiler/lib/intel64
	module unload compiler/gcc/4.8.5
	module unload gcc
	module load gcc-7.3.0/7.3.0
    fi
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${INTELROOT}
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${MKLROOT}
    export PATH=$PATH:${INTELROOT}
    export PATH=$PATH:${MKLROOT}

    module list

    NODES=$(cat $PBS_NODEFILE | cut -d " " -f 1 | sort | uniq | wc -l | cut -d " " -f 1)
fi
