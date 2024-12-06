import random
import cobra
from cobra.recon import *

model = cobra.io.read_sbml_model('/ibex/user/niuk0a/funcarve/cobra/bigg/iNJ661.xml') # 0.053
model = cobra.io.read_sbml_model('/ibex/user/niuk0a/funcarve/cobra/iYL1228.xml') # 1.043
print('model:1128:')
num_rm_rxn = int(float(0.05) * len(model.reactions))



rm_rxn = random.sample(model.reactions,num_rm_rxn)

model.remove_reactions(rm_rxn)
print('model1228:',len(model.reactions),model.optimize().objective_value,flush=True)
universe_model = cobra.io.load_json_model("/ibex/user/niuk0a/funcarve/cobra/bigg/universal_model_cobrapy.json")

pr2ec , predscore = read_clean_withscore('/ibex/user/niuk0a/funcarve/cobra/bigg/CP000647.1_sequence_iyj1228.csv')
universal_scoredict,r2maxecp = clean2biggr(predscore)
try:
    model.solver = 'glpk'
except:
    print('glpk not available')
try:
    model.solver = 'scip'
except:
    print('scip not available')
import time
# print('using penalties',flush=True)
# a = time.time()
# print('time:',time.time())
# solution = cobra.flux_analysis.gapfill(model,universe_model,penalties = universal_scoredict, demand_reactions=False)
# print('time:',time.time(),'time cost:',time.time()-a)
# print('reaction:',len(solution[0]))
# for reaction in solution[0]:
#     print(reaction.id)

# without penalties
print('without penalties',flush=True)
a = time.time()
print('time:',time.time())  
solution = cobra.flux_analysis.gapfill(model,universe_model, demand_reactions=False)
print('time:',time.time(),'time cost:',time.time()-a)
print('reaction:',len(solution[0]))
for reaction in solution[0]:
    print(reaction.id)
