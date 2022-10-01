from ase.clease import CEBulk , NewStructures , Concentration



# Database used to have all the structures
ECI_FILE = "eci_almgzn_hcp.json"
db_name = "almgzn_hcp_2x2x2.db"


def main():

    # Concentation arguments. NOTE: The way of specifying concentrations
    # will be changed in the future
##    conc_args = {
##        "conc_ratio_min_1": [[64, 0, 0]],
##        "conc_ratio_max_1": [[24, 40, 0]],
##        "conc_ratio_min_2": [[64, 0, 0]],
##        "conc_ratio_max_2": [[22, 21, 21]]
##    } # for a 4x4x4=64 FCC lattice
    
    conc_args = {
        "conc_ratio_min_1": [[8, 0, 0]],
        "conc_ratio_max_1": [[0, 8, 0]],
        "conc_ratio_min_2": [[6, 3, 4]],
        "conc_ratio_max_2": [[2, 5, 6]]
    } # for a 2x2x2=8 FCC lattice

##    conc_args = {
##        "conc_ratio_min_1": [[6, 6, 6]],
##        "conc_ratio_max_1": [[6, 6, 6]],
##        "conc_ratio_min_2": [[16, 6, 5]],
##        "conc_ratio_max_2": [[2, 12, 13]]
##    } # for a 3x3x2=18 HCP lattice


    basis_elements=[["Al" , "Mg" , "Zn"]]
    concentration = Concentration(basis_elements)
    ceBulk = CEBulk(
#        crystalstructure="hcp", a=3.21 , c=5.21 , size=[2, 2, 2],
        crystalstructure="hcp", a=3.21 , c=5.21 , size=[2, 2, 2],        
        concentration=concentration,
#        basis_elements=[["Al","Al", "Zn", "Zn"]], conc_args=conc_args,
        db_name=db_name, max_cluster_size=4
#        , max_cluster_dia=[0, 0, 5.3, 4.1, 4.1]
#        , max_cluster_dia=[0, 0, 5.3, 5.3, 5.3]
#        , max_cluster_dia=[0, 0, 7.0, 5.3, 5.3]
        , max_cluster_dia=[0, 0, 12.0, 4.5,5.3 ]
                                                    )    
    


    # Create an instance of the structure generator
    struc_generator = NewStructures(ceBulk, struct_per_gen=200)

    struc_generator.generate_initial_pool()


    struc_generator.generate_probe_structure(init_temp=1.0, final_temp=0.001,
    num_temp=5, num_steps_per_temp=100,approx_mean_var=True)
    
    ceBulk.view_clusters()
    # Generate new structures
##    struc_generator.generate_probe_structure()
##    ceBulk.view_clusters()

    # Evaluate and fit the ECIs
    # evaluate(ceBulk)

##def random_gs(bc, struct_gen):
##    from cemc.tools import GSFinder
##    with open(ECI_FILE, 'r') as infile:
##        eci = json.load(infile)
##
##    T = np.linspace(10, 2000, 30)[::-1]
##    composition = {"Al": 0.5, "Zn": 0.25, "Mg": 0.25}
##    gs_finder = GSFinder()
##    gs = gs_finder.get_gs(bc, eci, composition=composition, temps=T, num_steps_per_temp=1000)
##    print("Energy: {}".format(gs["energy"])
##    struct_gen.insert_structure(init_struct=gs["atoms"])
##


if __name__ == "__main__":
    main()
