from ase.clease import CEBulk as BulkCrystal
from ase.clease import Evaluate,  Concentration
from ase.clease import tools
import json


db_name = "almgzn_hcp_2x2x2_v2.db"


basis_elements=[["Al" , "Mg" , "Zn"]]
concentration = Concentration(basis_elements)
ceBulk = BulkCrystal(
crystalstructure="hcp", a=3.21, c=5.21, size=[2, 2, 2],
#crystalstructure="fcc", a=4.05, size=[3, 3, 3],
concentration=concentration,
db_name=db_name, max_cluster_size=4,
max_cluster_dia=[0, 0, 12.1, 8.6, 8.6]
)

tools.reconfigure(ceBulk,select_cond = [("converged","=",1)])

