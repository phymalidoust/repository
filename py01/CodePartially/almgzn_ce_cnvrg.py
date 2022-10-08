from ase.clease import CEBulk as BulkCrystal
from ase.clease import GenerateStructures
from ase.clease import Evaluate
import json

# File where ECIs are stored
eci_fname = "converged_v4/almgzn_hcp_2x2x2_MaxCD_4.json"



# Database used to have all the structures
db_name = "converged_v4/almgzn_hcp_2x2x2.db"


def main():

    # Concentation arguments. NOTE: The way of specifying concentrations
    # will be changed in the future
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
    
    
    ceBulk = BulkCrystal(
        crystalstructure="hcp", a=3.21 , c=5.21 , size=[2, 2, 2],
#        crystalstructure="fcc", a=4.05 , size=[2,2,2],        
        basis_elements=[["Al" , "Mg" , "Zn"]], conc_args=conc_args,
#        basis_elements=[["Al","Al", "Zn", "Zn"]], conc_args=conc_args,
        db_name=db_name, max_cluster_size=4
#        , max_cluster_dia=[0, 0, 5.3, 4.1, 4.1]
#        , max_cluster_dia=[0, 0, 5.3, 5.3, 5.3]
#        , max_cluster_dia=[0, 0, 7.0, 5.3, 5.3]
        , max_cluster_dia=[0, 0, 7.0, 6.5, 5.3]
                                                    )

    # Create an instance of the structure generator
#    struc_generator = GenerateStructures(ceBulk, struct_per_gen=20)

#    struc_generator.generate_initial_pool()
#    ceBulk.view_clusters()
    # Generate new structures
    #struc_generator.generate_probe_structure()

    # Evaluate and fit the ECIs
    evaluate(ceBulk)



def evaluate(BC):
    # Set up an Evaluator with L1 regularization
    evaluator = Evaluate(BC, fitting_scheme="l2")

    # Try different penalization value to find the best
    best_alpha = evaluator.plot_CV(1E-5, 1E-1, num_alpha=100)
    evaluator.set_fitting_scheme("l2" , best_alpha)

    # Find the ECIs using the best penalization value
    evaluator.plot_fit(interactive=False)
    print("Best penalization value: {}".format(best_alpha))
    eci_name = evaluator.get_cluster_name_eci(return_type="dict")

    with open(eci_fname, 'w') as outfile:
        json.dump(eci_name, outfile, indent=2, separators=(",", ":"))
    print("ECIs written to {}".format(eci_fname))


if __name__ == "__main__":
    main()
