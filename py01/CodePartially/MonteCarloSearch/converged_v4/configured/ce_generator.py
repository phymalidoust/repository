from ase.clease import CEBulk as BulkCrystal
from ase.clease import Evaluate, Concentration
import json
from ase.clease import ConvexHull
import matplotlib.pyplot as plt


set = 26

snglt = 4.5
dblt = 3.5
trplt = 3.5

start_alfa = 1E-5
end_alfa = 1E-2

num_alpha = 5000


# File where ECIs are stored
eci_fname = "set" + str(set) + "/" + "CEIs_fcc_MgZn.json"
fParaName = "set" + str(set) + "/" + "parameters.txt"






# Database used to have all the structures
db_name = "almgzn_hcp_2x2x2_v2.db"



basis_elements = [["Al", "Mg", "Zn"]]
concentration = Concentration(basis_elements)

ceBulk = BulkCrystal(
    crystalstructure="hcp", a=3.21, c=5.21, size=[2, 2, 2],
    # crystalstructure="fcc", a=4.05, size=[3, 3, 3],
        concentration=concentration,
        db_name=db_name, max_cluster_size=4,
        max_cluster_dia=[0, 0, snglt, dblt, trplt]
)



    
    # Set up an Evaluator with L1 regularization
evaluator = Evaluate(ceBulk, fitting_scheme="l1", select_cond=[("converged", "=", 1)])

    # Try different penalization value to find the best
best_alpha = evaluator.plot_CV(start_alfa, end_alfa, num_alpha=num_alpha)
evaluator.set_fitting_scheme("l1" , best_alpha)


ParamFile = open(fParaName, "w")
ParamFile.write("\nSingls diameter is: %f" % snglt)
ParamFile.write("\ndoubls diameter is: %f" % dblt)
ParamFile.write("\ntripels diameter is: %f" % trplt)
ParamFile.write("\nnumber of alpha steps is: %f" % num_alpha)


    # Find the ECIs using the best penalization value
evaluator.plot_fit(interactive=False)
print("Best penalization value: {}".format(best_alpha))
eci_name = evaluator.get_cluster_name_eci(return_type="dict")

with open(eci_fname, 'w') as outfile:
    json.dump(eci_name, outfile, indent=2, separators=(",", ":"))
print("ECIs written to {}".format(eci_fname))

chall =ConvexHull(ceBulk)
chall.plot()


