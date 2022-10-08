#!/usr/bin/env python3

from Model_Params import Model_Params
from Node import Node
import sys
import os

if __name__ == '__main__':
    nodes = 1
    try:
        nodes = int(sys.argv[1])
    except:
        pass

    cores = 32
    try:
        cores = int(sys.argv[2])
    except:
        pass

    ram = 122
    try:
        ram = int(sys.argv[3])
    except:
        pass

    sockets = 2
    try:
        sockets = int(sys.argv[4])
    except:
        pass

    max_cores_node = int(os.environ['BC_STANDARD_NODE_CORES'])
        
    node = Node(max_cores_node, ram, sockets)

    min_nkx = 16
    try:
        min_nkx = int(os.environ['spin_min_nkx'])
    except:
        pass
    
    twodims = False
    if (min_nkx == 0):
        twodims = True

    clustersize = 0

    # Find OMP_NUM_THREADS
    omp_threads = 1
    if(node.cores % 7 == 0):
        omp_threads = 7
    else:
        for i in range(5,1,-1):
            if(node.cores % i == 0):
                omp_threads=i
                break

    # total number of ranks
    total_ranks = cores
    ranks_per_node = int(cores / nodes)       

    if (ranks_per_node * omp_threads > node.cores):
        omp_threads = int(node.cores/ranks_per_node)

    # number of rows and cols in each processor grid
    pg_rows = 1
    done = False
    for j in range(6,0,-1):
        if done:
            break
        if(total_ranks % j == 0):
            pg_rows = j
        else:
            continue

        model = Model_Params()
        model.read_env()

        if (twodims):
            range_start = int(total_ranks/pg_rows)
            range_end = range_start - 1
            range_incr = -1
            if (range_start < range_end):
                continue
        else:
            range_start = pg_rows
            range_end = int(total_ranks/pg_rows)+1
            range_incr = 1

        pg_cols = 1

        for i in range(range_start, range_end, range_incr):
            grid_size = i*pg_rows
            if(total_ranks % grid_size != 0):
                continue
    
            check = model.size_check(pg_rows, i, clustersize, twodims)
            if(not check[0]):
                continue

            if(model.ram_per_rank * ranks_per_node > 0.85 * ram):
                continue

            if(check[1] and check[2]):
                pg_cols = i
                done = True
                break

            while(False in check):
                model.kx_factor += 0.001

                check = model.size_check(pg_rows, i, clustersize, twodims)

                if(not check[0]):
                    break

                if((not twodims) and (model.kx_factor > 2.0)):
                    break

                if ((twodims) and (model.kx_factor > 3.0)):
                    break

                if(model.ram_per_rank * ranks_per_node > 0.95 * ram):
                    break
   
            if(     not (False in check) and \
                (model.ram_per_rank * ranks_per_node < 0.95 * ram)):
                pg_cols = i
                done = True
                break
            else:
                model.kx_factor = model.min_kx_factor

    check = model.size_check(pg_rows,pg_cols,clustersize, twodims)
 
    if (not twodims) and (False in check):
        # Despite our best efforts we couldn't create an optimum arrangment
        # Reset PG to 1x1 and allow spawning of additional nkx
        model.n_proc_rows = 1
        model.n_proc_cols = 1
    
    if (twodims) and (pg_rows*pg_cols != total_ranks):
        # This isn't going to work and we cannot allow more nkx in 2D cases
        # The only way to signal failure out of this script is to assign
        # an illegal value to omp_threads and have the PBS script test for it
        omp_threads = 100
        
    grid_size = model.n_proc_cols * model.n_proc_rows
    flight = int(total_ranks / grid_size)
    n_flights = int((min_nkx + flight - 1) / flight)
    model.nkx = n_flights * flight
    # If we have set min_nkx to 0 to signal we desire a 2D run,
    # model.nkx will now still be zero. That won't work, we need
    # to override to 1 so the code will run
    if twodims:
        model.nkx = 1

    print(f'total_ranks : {total_ranks}')
    print(f'ranks/node  : {ranks_per_node}')
    print(f'n_proc_rows : {model.n_proc_rows}')
    print(f'n_proc_cols : {model.n_proc_cols}')
    print(f'flight_size : {flight}')
    print(f'# of flights: {n_flights}')
    print(f'nkx         : {model.nkx}')
    print(f'Ram per rank: {model.ram_per_rank:.2f} Gb')
    print(f'kx_factor   : {model.kx_factor:.3f}')
    print(f'omp_threads : {omp_threads}')

    if (omp_threads == 100):
        print("\n# of cores isn't compatible with 2D problem conditions")
        exit()

    model.write_input('input.txt')
