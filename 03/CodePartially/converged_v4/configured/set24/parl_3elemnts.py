from __future__ import print_function
from ase.clease import CEBulk
from ase.clease import Evaluate
from ase.clease import Concentration
from cemc import get_atoms_with_ce_calc
import json
import numpy as np

import os
import sys

snglt = 4.5
dblt = 4.5
trplt = 3.5

# Database used to have all the structures
ECI_FILE = "CEIs_fcc_MgZn.json"
db_name = "almgzn_hcp_2x2x2_v2.db"

my_rank = int(sys.argv[1])
num_ranks = int(sys.argv[2])

##if my_rank==0:
##    if not os.path.exists('data/energies'):
##        os.makedirs('data/energies')
##    if not os.path.exists('data/largstrctrs'):
##        os.makedirs('data/largstrctrs')

##folder1 = "/data/energies"
##folder2 = "/data/largstrctrs"

folder1 = "./"
folder2 = "./"

pname = str(my_rank)
fname = folder1 + '/' + 'cpu_' + pname + '.txt'
Lg_dbname = folder2 + '/' + 'cpu_' + pname + '_large_10x10x10.db'
tname = folder1 + '/' + 'Alvc.txt'


# fname = 'file.txt'
# Lg_dbname = 'larg.db'


def main():
    # Concentation arguments. NOTE: The way of specifying concentrations
    # will be changed in the future
    basis_elements = [["Al", "Mg", "Zn"]]
    conc = Concentration(basis_elements)

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
        # "crystalstructure": "fcc", "a": 4.05, "size": [3, 3, 3],
              "crystalstructure":"hcp", "a":3.21, "c":5.21, "size":[2, 2, 2],
        "concentration": conc,
        "db_name": db_name, "max_cluster_size": 4
        #        , "max_cluster_dia":[0, 0, 5.3, 4.1, 4.1]
        #        , "max_cluster_dia":[0, 0, 5.3, 5.3, 5.3]
        , "max_cluster_dia": [0, 0, snglt, dblt, trplt]
        #        , "max_cluster_dia":[0, 0, 7.0, 6.5, 5.3]
    }

    #    wait = input("PRESS ENTER TO CONTINUE.")

    ceBulk = CEBulk(**argms)

    #    wait = input("PRESS ENTER TO CONTINUE.")

    ceBulk.reconfigure_settings()
    # ceBulk.reconfig_db_entries()

    mc_cell_size = [10, 10, 10]

    with open(ECI_FILE, 'r') as infile:
        eci = json.load(infile)

    ##    eci = {"c1_0": 0.1}

    atoms = get_atoms_with_ce_calc(ceBulk, argms, eci=eci, size=mc_cell_size, db_name=Lg_dbname)

    calc = atoms.get_calculator()
    ceBulk = calc.BC
    ceBulk.atoms.set_calculator(calc)

    # Create an instance of the structure generator
    # struc_generator = GenerateStructures(ceBulk, struct_per_gen=4)

    # struc_generator.generate_initial_pool()
    # ceBulk.view_clusters()
    # Generate new structures
    # struc_generator.generate_probe_structure()

    # Evaluate and fit the ECIs
    # evaluate(ceBulk)

    ##				comm.Barrier()

    alvc = np.arange(0.05, 0.98, 0.0207)

    strt = my_rank * len(alvc) / num_ranks
    endd = (my_rank + 1) * len(alvc) / num_ranks

    # strt =  1
    # endd = 10

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
        znvc = np.arange(0.05, 1.0 - al0 - 0.01, 0.0207)
        for nzn in range(0 , len(znvc)):
            zn0 = znvc[nzn]
            mg0 = 1 - al0 - zn0
            enrg = random_gs(ceBulk, al0, zn0, mg0)
            fd = open(fname , 'a+')
            print(enrg, al0, zn0, mg0, file=fd)
            fd.close()


def random_gs(ceBulk, aal, zzn, mmg):
    from cemc.tools import GSFinder
    with open(ECI_FILE, 'r') as infile:
        eci = json.load(infile)

    ##    eci = {"c1_0": 0.1}

    T = np.linspace(10, 2000, 700)[::-1]
    composition = {"Al": aal, "Zn": zzn, "Mg": mmg}
    gs_finder = GSFinder()
    gs = gs_finder.get_gs(ceBulk, eci, temps=T, n_steps_per_temp=2500, composition=composition)
    # gs = gs_search.get_gs(bc, eci, temps=temps, n_steps_per_temp=100, composition=composition)
    # print("Energy: {}".format(gs["energy"]))
    ##    try:
    ##        struct_gen.insert_structure(init_struct=gs["atoms"])
    ##    except RuntimeError as exc:
    ##        print(str(exc))
    ##        pass

    return gs["energy"]


if __name__ == "__main__":
    main()
