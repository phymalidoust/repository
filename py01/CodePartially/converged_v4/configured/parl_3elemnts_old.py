from __future__ import print_function
from ase.clease import CEBulk
from ase.clease import GenerateStructures
from ase.clease import Evaluate
from cemc import get_ce_calc
import json
import numpy as np

import os
import sys


# Database used to have all the structures
ECI_FILE = "almgzn_hcp_2x2x2.json"
db_name = "almgzn_hcp_2x2x2.db"


my_rank   = int(sys.argv[1])
num_ranks = int(sys.argv[2])

##if my_rank==0:
##    if not os.path.exists('data/energies'):
##        os.makedirs('data/energies')
##    if not os.path.exists('data/largstrctrs'):
##        os.makedirs('data/largstrctrs')        

folder1 = "./"
folder2 = "./"
pname = str(my_rank)
fname = folder1 + '/' + 'cpu_' + pname + '.txt'
Lg_dbname = folder2 + '/' + 'cpu_' + pname + '_large_4x4x4.db'
tname  = folder1 + '/' + 'Alvc.txt'

def main():

    # Concentation arguments. NOTE: The way of specifying concentrations
    # will be changed in the future
    conc_args = {
        "conc_ratio_min_1": [[8, 0, 0]],
        "conc_ratio_max_1": [[0, 8, 0]],
        "conc_ratio_min_2": [[6, 3, 4]],
        "conc_ratio_max_2": [[2, 5, 6]]
    }
    
##    argms = {
##        "crystalstructure":"fcc", "a":4.05, "size":[10, 10, 10],
##        "basis_elements":[["Al", "Mg", "Zn"]], "conc_args":conc_args,
##        "db_name":db_name
##    }


##    argms = {
##        "crystalstructure":"fcc", "a":4.05, "size":[4, 4, 4],
##        "basis_elements":[["Al", "Mg", "Zn"]], "conc_args":conc_args,
##        "db_name":db_name, "max_cluster_size":4
##        , "max_cluster_dia":[0, 0, 7.0, 7.0, 6.2]
##    }

##			for r in range(0,num_ranks-1):
##			
##					if r==my_rank:
    argms = {
 #       "crystalstructure":"fcc", "a":4.05, "size":[4, 4, 4],
        "crystalstructure":"hcp", "a":3.21, "c":5.21, "size":[2, 2, 2],
        "basis_elements":[["Al", "Mg", "Zn"]], "conc_args":conc_args,
        "db_name":db_name, "max_cluster_size":4
#        , "max_cluster_dia":[0, 0, 5.2, 4.1, 4.1]
        , "max_cluster_dia":[0, 0, 7.0, 6.5, 5.3]
    }



    ceBulk = CEBulk( **argms)

    #ceBulk.reconfigure_settings()
    #ceBulk.reconfig_db_entries()

    mc_cell_size = [10, 10, 10]

    with open(ECI_FILE, 'r') as infile:
        eci = json.load(infile)

##    eci = {"c1_0": 0.1}

    calc = get_ce_calc(ceBulk, argms, eci=eci, size=mc_cell_size, db_name=Lg_dbname)

    ceBulk = calc.BC
    ceBulk.atoms.set_calculator( calc )
    
    

    
    # Create an instance of the structure generator
    struc_generator = GenerateStructures(ceBulk, struct_per_gen=4)

    #struc_generator.generate_initial_pool()
    #ceBulk.view_clusters()
    # Generate new structures
    #struc_generator.generate_probe_structure()
    
    # Evaluate and fit the ECIs
    # evaluate(ceBulk)

##				comm.Barrier()

    alvc = np.arange(0.05, 0.95, 0.019) # makes 48 steps. Then if we request 8 processors, each will do 6 steps.
                                        # if request 24 processors, each will do 2 steps
    
    strt =  my_rank*len(alvc)/num_ranks
    endd = (my_rank + 1)*len(alvc)/num_ranks
        
##    if my_rank == 0:
##        strt =  my_rank*len(alvc)/num_ranks 
        
##    print('start: ', strt, 'endd: ', endd )
##    print('Range: ', range(int(strt) , int(endd)))
##    wait = input("PRESS ENTER TO CONTINUE.")

##    for nal in range(0,len(alvc)):
##        al0 = alvc[nal]
##        znvc = np.arange(0.03, 1.0 - al0 - 0.01, 0.025)
##        for nzn in range(0 , len(znvc)):
##            zn0 = znvc[nzn]
##            mg0 = 1 - al0 - zn0 
##            fd = open(tname , 'a+')
##            print(al0,zn0,mg0,file=fd)
##            fd.close()

    for nal in range(int(strt) , int(endd)):
        al0 = alvc[nal]
        znvc = np.arange(0.03, 1.0 - al0 - 0.01, 0.025)
        for nzn in range(0 , len(znvc)):
            zn0 = znvc[nzn]
            mg0 = 1 - al0 - zn0 
            enrg = random_gs(ceBulk , struc_generator , al0 , zn0 , mg0)
            fd = open(fname , 'a+')
            print(enrg,al0,zn0,mg0,file=fd)
            fd.close()

def random_gs(bc, struct_gen , aal , zzn , mmg):
    from cemc.tools import GSFinder
    with open(ECI_FILE, 'r') as infile:
        eci = json.load(infile)

##    eci = {"c1_0": 0.1}

    T = np.linspace(200, 2000, 500)[::-1]
    composition = {"Al": aal, "Zn": zzn, "Mg": mmg}
    gs_finder = GSFinder()
    gs = gs_finder.get_gs(bc, eci, composition=composition, temps=T, n_steps_per_temp=2500)
    #print("Energy: {}".format(gs["energy"]))
##    try:
##        struct_gen.insert_structure(init_struct=gs["atoms"])
##    except RuntimeError as exc:
##        print(str(exc))
##        pass

    return gs["energy"]



if __name__ == "__main__":
    main()
