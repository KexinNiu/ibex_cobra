######## model iAF987 gapfilling ########

from copy import deepcopy
import warnings
import cobra
from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis import gapfill
from cobra.io import load_model
from utils import *
import pandas as pd 
import random
from cobra.recon import _checkModel, _find_reactions, _gapfill_model, _add_annotation, w_find_reactions,print_model_info,read_clean_withscore,clean2biggr,_set_base_inputs
import mergem
from mergem import translate, load_model, save_model
import argparse
import re
parser = argparse.ArgumentParser(description='Precentage to remove reactions')
parser.add_argument('--prec', default = 0.05, help='percentage of reactions to remove')
parser.add_argument('--seednum', default = 1234, help='seed for random')
parser.add_argument('--type', default = 'all', help='seed for random')
args = parser.parse_args()
print('args:',args.prec,args.seednum,args.type)
result_log = {}
result_log['prec'] = args.prec
result_log['seednum'] = args.seednum
result_log['type'] = args.type
model = cobra.io.read_sbml_model('/ibex/user/niuk0a/funcarve/cobra/iAF987.xml')
solution = model.optimize()
print_model_info(model) 
result_log['originalFBA'] = solution.objective_value
print('model.media:',model.medium)

names,seqs = read_fasta("/ibex/user/niuk0a/funcarve/cobra/bigg/aaseq_CP000148.1.txt")
pr2ec , predscore = read_clean_withscore('/ibex/user/niuk0a/funcarve/cobra/aaseq_CP000148.1_maxsep.csv')

unreviewed_f = '/ibex/user/niuk0a/funcarve/cobra/269799_unreviewed.tsv'
unreviewed = pd.read_csv(unreviewed_f, sep='\t')
print(len(unreviewed))
unreviewed_genes = set(unreviewed['Gene Names'].str.split(' ').sum())
print(len(unreviewed_genes))

def get_rm_pool(model,unreviewed_genes,type='all'):
    can_remove = []
    genes_rxn = {}
    rxn_genes = {}
    for rxn in model.reactions:
        if rxn.gene_reaction_rule != '':
            if 'or' in rxn.gene_reaction_rule or 'and' in rxn.gene_reaction_rule:
                if 'and' in rxn.gene_reaction_rule:
                    genes = rxn.gene_reaction_rule.split(' and ')
                    # genes_rxn[genes] = rxn
                    for gene in genes:
                        if gene in genes_rxn:
                            genes_rxn[gene].append(rxn.id)
                        else:
                            genes_rxn[gene] = [rxn.id]
                    
                        if  gene in unreviewed_genes:
                            can_remove.append(rxn)
                            break
                    rxn_genes[rxn.id] = "|".join(genes)
                elif type == 'one':
                    genes = rxn.gene_reaction_rule.split(' or ')
                    rxn_genes[rxn.id] = genes
                    for gene in genes:
                        # rxn_genes[rxn.id] = genes
                        if gene in genes_rxn:
                            genes_rxn[gene].append(rxn.id)
                        else:
                            genes_rxn[gene] = [rxn.id]
                        if gene in unreviewed_genes:
                            can_remove.append(rxn)
                            break
                elif type == 'all':
                    genes = rxn.gene_reaction_rule.split(' or ')
                    rxn_genes[rxn.id] = genes
                    flag = 0
                    for gene in genes:
                        # rxn_genes[rxn.id] = genes
                        if gene in genes_rxn:
                            genes_rxn[gene].append(rxn.id)
                        else:
                            genes_rxn[gene] = [rxn.id]
                        
                        if gene not in unreviewed_genes:
                            flag = 1
                    if flag == 0:
                        can_remove.append(rxn)
    return can_remove,genes_rxn,rxn_genes


can_remove,genes_rxn,rxn_genes = get_rm_pool(model,unreviewed_genes,type=args.type)



num_rm_rxn = int(float(args.prec) * len(can_remove))
random.seed(int(args.seednum))
rm_rxn = random.sample(can_remove, num_rm_rxn)
rm_rxn_id = [r.id for r in rm_rxn]
print(f"total reactions: {len(model.reactions)}")
print(f"can remove reactions: {len(can_remove)}")
print(f"5% rm_rxn: {len(rm_rxn)}")
result_log['ori_model_reaction'] = len(model.reactions)
result_log['ori_model_metabolites'] = len(model.metabolites)
result_log['ori_model_genes'] = len(model.genes)
result_log['rm_reaction'] = len(rm_rxn)

universal = cobra.io.load_json_model("/ibex/user/niuk0a/funcarve/cobra/bigg/universal_model_cobrapy.json")
print('len(universal.reactions)',len(universal.reactions))
# universal_scoredict = clean2biggr(universal,predscore)
universal_scoredict = clean2biggr(predscore)

model.remove_reactions(rm_rxn)
model,rm_metabolits = cobra.manipulation.prune_unused_metabolites(model)

result_log['rm_metabolites'] = len(rm_metabolits)
result_log['afrm_metabolites'] = len(model.metabolites)
result_log['afrm_reactions'] = len(model.reactions)
result_log['uni_reactions'] = len(universal.reactions)

for r in rm_rxn:
    if r not in universal.reactions:
        print('r:',r.id,'not in universal')
        universal.add_reactions([r.copy()])

exchange_arg = 1 # 0 for no exchange, 1 for exchange
metabolic_tasks=[] # no metabolic tasks
min_frac = 0.01 # Minimum objective fraction required during gapfilling
max_frac = 0.5  # Maximum objective fraction required during gapfilling
file_type = 1 # 3 for sbml draft model (remove reactions from model can not be ==3
gram_type ='negative'
media = model.medium
print_model_info(model)
gfmodel = deepcopy(model)
# for gfmodel obj translation

universal_obj = str(gfmodel.objective.expression).split()[0].split('*')[-1]
# reac_univ_id = mergem.map_reaction_univ_id(universal_obj)
# transobj = mergem.get_reaction_properties(reac_univ_id)['ids']
# flag = 0
# for item in transobj:
#     if 'seed' in item:
#         flag = 1
#         gf_universal_obj = item.split(':')[1]
#         gfmodel.objective = gfmodel.reactions.get_by_id(universal_obj)
#         break
# if flag == 0:
#     print('no seed objective')
#     gf_universal_obj = 'biomass_GmNeg'
# # end of obj translation   

gfmodel = mergem.translate(gfmodel,'seed')

gfwmodel = deepcopy(model)  
# biomass_re = re.compile("biomass", re.IGNORECASE)
# possible_objectives = m.reactions.query(biomass_re)


################### gapfilling NO weighted #####################
import pickle
## bigg database
from reframed import load_cbmodel, save_cbmodel
# universe_file ='/ibex/user/niuk0a/funcarve/carveme/data/generated/universe_bacteria.xml.gz'
# universe_model = load_cbmodel(universe_file, flavor='bigg')
universe_model = cobra.io.load_json_model("/ibex/user/niuk0a/funcarve/cobra/bigg/universal_model_cobrapy.json")
pr2ec , predscore = read_clean_withscore('/ibex/user/niuk0a/funcarve/cobra/aaseq_CP000148.1_maxsep.csv')
universal_scoredict,r2maxecp = clean2biggr(predscore)
model = cobra.io.read_sbml_model('/ibex/user/niuk0a/funcarve/cobra/iAF987_555_0.05_rm.xml')
solution = cobra.flux_analysis.gapfill(model, universe_model,penalties = universal_scoredict, demand_reactions=False)
for reaction in solution[0]:
    print(reaction.id)
'''
ModelSEED 
import pickle
universal_seedf = '/ibex/user/niuk0a/funcarve/cobra/universal_reconori.pickle'
universal_seed = pickle.load(open(universal_seedf,'rb'))
print('################### gapfilling NO weighted #####################')
# universal_obj =
print('Identifying new metabolism...')
draft_reactions = set([x.id for x in gfmodel.reactions])
draft_metabolites = set([x.id for x in gfmodel.metabolites])
warnings.filterwarnings('ignore')
new_reactions = _find_reactions(gfmodel, universal_seed, metabolic_tasks, gf_universal_obj, min_frac, max_frac, 1, file_type)
print(new_reactions)
inter = new_reactions.intersection(set(rm_rxn_id))
result_log['addgf_media_reactions_inter'] = len(inter)  
print('gf_newR:',len(new_reactions),'intersection:',len(inter))
filled_genre = _gapfill_model(gfmodel, universal_seed, new_reactions, gf_universal_obj, 1)
if file_type != 3:
    print('Identifying new metabolism (Step 2 of 2)...')
    filled_genre = _set_base_inputs(filled_genre, universal_seed)
    media_reactions = _find_reactions(filled_genre, universal_seed, metabolic_tasks, gf_universal_obj, min_frac, max_frac, 2, file_type)
    final_genre = _gapfill_model(filled_genre, universal_seed, media_reactions, universal_obj, 2)
    final_genre = _add_annotation(final_genre, gram_type)
else: 
    final_genre = _add_annotation(filled_genre, gf_universal_obj)
final_genre = _add_annotation(filled_genre, gf_universal_obj)
print('GF new reactions:',len(new_reactions))
print('GF media reactions:',len(media_reactions))
allnewid = set(new_reactions).union(set(media_reactions))
inter = allnewid.intersection(set(rm_rxn_id))
result_log['addgf_media_reactions_inter'] = len(inter)  
print('gf_newR:',len(new_reactions),'intersection:',len(inter))
# Correct exchanges and check new model
if exchange_arg == 0:
    for exch in final_genre.exchanges: exch.bounds = (0., 0.)
else:
    for exch in final_genre.exchanges: exch.bounds = (-1000., 1000.)
for rxn in final_genre.reactions:
    if 'Exchange reaction for' in rxn.name:
        rxn.name = list(rxn.metabolites)[0].name + ' exchange'
biomass = _checkModel(draft_reactions, draft_metabolites, final_genre)
print_model_info(final_genre)
sol = final_genre.optimize()
result_log['gf_FBA'] = sol.objective_value
result_log['gf_FBA_status'] = sol.status
result_log['gf_metabolites'] = len(final_genre.metabolites)
result_log['gf_reactions'] = len(final_genre.reactions)
result_log['gf_genes'] = len(final_genre.genes)
## save the final model
cobra.io.write_sbml_model(final_genre, f'/ibex/user/niuk0a/funcarve/cobra/iAF987_{str(args.seednum)}_{str(args.prec)}_{args.type}_gf.xml')
print('finish save gf model',flush=True)'''







################### gapfilling weighted #####################
universal = cobra.io.load_json_model("/ibex/user/niuk0a/funcarve/cobra/bigg/universal_model_cobrapy.json")
if len(media) != 0:
    media_condition = set(['EX_' + cpd for cpd in media])
    universal_reactions = set([x.id for x in universal.reactions])
    for rxn in universal_reactions:
        if rxn.startswith('EX_') == True:
            universal.reactions.get_by_id(rxn).bounds = (0, 1000.0)
        if rxn in media_condition:
            universal.reactions.get_by_id(rxn).bounds = (-1000.0, 10000)
# universal_obj = str(model.objective.expression).split()[0].split('*')[-1]
print('################### gapfilling weighted #####################')
universal_obj = str(gfwmodel.objective.expression).split()[0].split('*')[-1]
pr2ec , predscore = read_clean_withscore('/ibex/user/niuk0a/funcarve/cobra/aaseq_CP000148.1_maxsep.csv')
universal_scoredict,r2maxecp = clean2biggr(predscore) # r2maxecp[rid]=[pr,ec,score]
print('Identifying new metabolism...')
draft_reactions = set([x.id for x in model.reactions])
draft_metabolites = set([x.id for x in model.metabolites])
warnings.filterwarnings('ignore')
new_reactions_ids,coefficientDict = w_find_reactions(model, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 1, file_type,universal_scoredict)
if file_type != 3:
    print('Identifying new metabolism (Step 2 of 2)...')
    filled_genre = _set_base_inputs(filled_genre, universal)
    media_reactions,coefficientDict = w_find_reactions(filled_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 2, file_type,universal_scoredict)
    final_genre = _gapfill_model(filled_genre, universal, media_reactions, universal_obj, 2)
    final_genre = _add_annotation(final_genre, gram_type)
else: 
    final_genre = _add_annotation(filled_genre, universal_obj)
print('GFW new reactions:',len(new_reactions_ids))
print('GFW media reactions:',len(media_reactions))
allnewid = set(new_reactions_ids).union(set(media_reactions))
inter = allnewid.intersection(set(rm_rxn_id))
result_log['addgfw_reactions_inter'] = len(inter)  
print('gf_newR:',len(new_reactions_ids),'intersection:',len(inter))
for rid in new_reactions_ids:
    print(f"reaction:{rid} add with {r2maxecp[rid]}")


filled_genre = _gapfill_model(model, universal, new_reactions_ids, universal_obj, 1)
final_genre = _add_annotation(filled_genre, universal_obj)

# Correct exchanges and check new model
if exchange_arg == 0:
    for exch in final_genre.exchanges: exch.bounds = (0., 0.)
else:
    for exch in final_genre.exchanges: exch.bounds = (-1000., 1000.)
for rxn in final_genre.reactions:
    if 'Exchange reaction for' in rxn.name:
        rxn.name = list(rxn.metabolites)[0].name + ' exchange'
biomass = _checkModel(draft_reactions, draft_metabolites, final_genre)
sol = final_genre.optimize()
blockr = cobra.flux_analysis.find_blocked_reactions(final_genre)
result_log['Wafgf_rm_blocked_reactions'] = len(blockr)
result_log['W_FBA'] = sol.objective_value
result_log['W_FBA_status'] = sol.status
result_log['W_metabolites'] = len(final_genre.metabolites)
result_log['W_reactions'] = len(final_genre.reactions)
result_log['W_genes'] = len(final_genre.genes)
print_model_info(final_genre)

## save the final model
cobra.io.write_sbml_model(final_genre, f'/ibex/user/niuk0a/funcarve/cobra/iAF987_{str(args.seednum)}_{str(args.prec)}_{args.type}_gfw.xml')
print('finish save gfw model',flush=True)



outf = f'result_log_iAF987_{args.seednum}_allunr_orirecon.txt'
print('outf:',outf) 
with open(outf,'a') as ouf:
    print('#######################',file=ouf)
    print(f'# precetage = {args.prec} #',file=ouf)
    print('\t'.join([str(k) for k in result_log.keys()]),file=ouf)
    print('\t'.join([str(v) for v in result_log.values()]),file=ouf)
    print('#######################',file=ouf)
ouf.close()
print('#######################',flush=True)
print(f'# precetage = {args.prec} #',flush=True)
print('\t'.join([str(k) for k in result_log.keys()]),flush=True)
print('\t'.join([str(v) for v in result_log.values()]),flush=True)
print('#######################',flush=True)
print(f'done with {args.prec}',flush=True)

''''

Traceback (most recent call last):
  File "/ibex/user/niuk0a/funcarve/cobra/testgapfill.py", line 171, in <module>
    filled_genre = _gapfill_model(gfmodel, universal_seed, new_reactions, gf_universal_obj, 1)
  File "/ibex/user/niuk0a/funcarve/cobra/recon.py", line 127, in _gapfill_model
    model.add_boundary(cpd, type='exchange', reaction_id=exch_id, lb=-1000.0, ub=1000.0)
  File "/ibex/user/niuk0a/anaconda3/envs/cobra/lib/python3.9/site-packages/cobra/core/model.py", line 609, in add_boundary
    raise ValueError(
ValueError: The metabolite is not an external metabolite (compartment is `extracellular` but should be `e`). Did you mean to add a demand or sink? If not, either change its compartment or rename the model compartments to fix this.
''''