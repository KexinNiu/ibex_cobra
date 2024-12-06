import cobra
from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis import gapfill
from cobra.io import load_model
from utils import *
import pandas as pd 
import random
from recon_gf import weighted_find_reactions, _add_annotation


import argparse
parser = argparse.ArgumentParser(description='Precentage to remove reactions')
parser.add_argument('--prec', default = 0.05, help='percentage of reactions to remove')
parser.add_argument('--seednum', default = 1234, help='seed for random')
parser.add_argument('--model', default = '/ibex/user/niuk0a/funcarve/cobra/iAF987.xml', help='model file')
parser.add_argument('--genome', default = '/ibex/user/niuk0a/funcarve/cobra/bigg/aaseq_CP000148.1.txt', help='genome file')
parser.add_argument('--clean', default = '/ibex/user/niuk0a/funcarve/cobra/aaseq_CP000148.1_maxsep.csv', help='clean file')
parser.add_argument('--unreviewed', default = '/ibex/user/niuk0a/funcarve/cobra/269799_unreviewed.tsv', help='unreviewed file')
parser.add_argument('--gf', default = 'w', help='unreviewed file')

args = parser.parse_args()
print('args:',args.prec,args.seednum)
result_log = {}
#### load model ####
modelf = '/ibex/user/niuk0a/funcarve/cobra/iAF987.xml'
print('model:',modelf)
model = cobra.io.read_sbml_model(modelf)
solution = model.optimize()
print('#ori model FBA solution:',solution)
result_log['ori_model FBA'] = solution.objective_value
print(model.medium)

### load genome data ##
## matching genome CP000148.1 
genomef = "/ibex/user/niuk0a/funcarve/cobra/bigg/aaseq_CP000148.1.txt"
print('genomef:',genomef)
names,seqs = read_fasta(genomef)
cleanf = '/ibex/user/niuk0a/funcarve/cobra/aaseq_CP000148.1_maxsep.csv'
print('cleanf:',cleanf)

pr2ec , predscore = read_clean_withscore(cleanf)
with open('biggr2ec.pkl', 'rb') as f:
    biggr2ec = pickle.load(f)
with open('biggec2r.pkl', 'rb') as f:
    biggec2r = pickle.load(f)

## load unreviewed proteins
unreviewed_f = '/ibex/user/niuk0a/funcarve/cobra/269799_unreviewed.tsv'
print ('unreviewed_f:',unreviewed_f)
unreviewed = pd.read_csv(unreviewed_f, sep='\t')

unreviewed_genes = set(unreviewed['Gene Names'].str.split(' ').sum())
print('unreviewed',len(unreviewed),'unreviewed_genes',len(unreviewed_genes))


## select reactions to remove ## 
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

            # else:
            #     genes = rxn.gene_reaction_rule.split(' or ')
            #     rxn_genes[rxn.id] = genes
            #     for gene in genes:
            #         # rxn_genes[rxn.id] = genes
            #         if gene in genes_rxn:
            #             genes_rxn[gene].append(rxn.id)
            #         else:
            #             genes_rxn[gene] = [rxn.id]
            #         if gene in unreviewed_genes:
            #             can_remove.append(rxn)
            #             break
            else:
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

# random remove 5% reactions in can_remove
# num_rm_rxn = int(len(can_remove) * 0.05)
num_rm_rxn = int(float(args.prec) * len(can_remove))
# seed number
random.seed(int(args.seednum))

rm_rxn = random.sample(can_remove, num_rm_rxn)
rm_rxn_id = [r.id for r in rm_rxn]
print(f"total reactions: {len(model.reactions)}")
print(f"can remove reactions: {len(can_remove)}")
print(f"5% rm_rxn: {len(rm_rxn)}")
# exit()
result_log['ori_model reaction'] = len(model.reactions)
result_log['ori_model metabolites'] = len(model.metabolites)
result_log['ori_model genes'] = len(model.genes)
result_log['rm_reaction'] = len(rm_rxn)



# universal = cobra.io.load_json_model("/ibex/user/niuk0a/funcarve/cobra/bigg/universal_model_cobrapy.json")
universal = cobra.io.read_sbml_model("/ibex/user/niuk0a/funcarve/cobra/universal_model_cobrapy_filter.xml")

print('len(universal.reactions)',len(universal.reactions))
# universal_scoredict = clean2biggr(universal,predscore)
universal_scoredict = clean2biggr(predscore)

print('rm_rxn in r2ec',len([r for r in rm_rxn if r.id in biggr2ec.keys()]),'total:',len(rm_rxn))

# remove reactions from model
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

oriobj = model.objective
for r in model.reactions:
    if 'BIOMASS' in r.id:
        print('r|',r.id,sep='')
        model.objective = r
if oriobj == model.objective:
    print('!!no biomass reaction found')
else:
    print('biomass reaction found:',model.objective)

print('model.objective.expression:',model.objective.expression)
model.solver.configuration.tolerances.feasibility = 1e-6
# model.solver.configuration.tolerances.optimality = 1e-6

media = model.medium.keys()
######## media ########
######################################################
########### weighted pFBA gapfiller from recon_gf.py###########
##############################################
from recon_gf import _find_reactions, _gapfill_model, _set_base_inputs
draft_genre = model
# save draft model
# cobra.io.write_sbml_model(draft_genre, f'iAF987_{args.seednum}_{args.prec}_rm.xml')

metabolic_tasks = [] # not required
# universal_obj = model.objective
universal_obj = str(draft_genre.objective.expression).split()[0].split('*')[-1]
print('universal_obj:',universal_obj)   

min_frac = 0.01
max_frac = 0.5
file_type = 1
if args.gf == 'w':

    universal.add_reactions([model.reactions.get_by_id(universal_obj)])
    new_reactions ,coefficientDict = weighted_find_reactions(draft_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 1, file_type, universal_scoredict,rm_rxn)
    print('@@@ this is weighted new_reactions ids:',new_reactions)
    result_log['addW_reactions'] = len(new_reactions)
    # reaction are removed from model
    inter = set(rm_rxn_id).intersection(set(new_reactions)) 
    print('intersection:',inter,'total added:',len(new_reactions),'removed reactions:',len(rm_rxn))
    result_log['addW_reactions_inRM'] = inter
    result_log['addW_reactions_inRM_num'] = len(inter)

    filled_genre = _gapfill_model(draft_genre, universal, new_reactions, universal_obj, 1)
    print('filled_genre:',len(filled_genre.reactions))
    sol = filled_genre.optimize()
    print('fba weighted gapped filled_genre model:',sol.objective_value)

    if file_type != 3:
        print('Identifying new metabolism (Step 2 of 2)...')
        filled_genre = _set_base_inputs(filled_genre, universal)
        media_reactions,coefficientDict = weighted_find_reactions(filled_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 2, file_type,universal_scoredict,rm_rxn)
        print('media_reactions:',len(media_reactions))
        final_genre = _gapfill_model(filled_genre, universal, media_reactions, universal_obj, 2)
        print('final_genre:',len(final_genre.reactions))
    else:
        final_genre = _add_annotation(filled_genre, universal_obj)
    alladd = set(new_reactions).union(media_reactions)
    inter = alladd.intersection(set(rm_rxn_id))
    print('alladd:',len(alladd),'intersection:',len(inter))
    result_log['addgf_media_reactions_inter'] = len(inter)    

    sol = final_genre.optimize()
    blockr = cobra.flux_analysis.find_blocked_reactions(final_genre)
    result_log['Wafgf_rm_blocked_reactions'] = len(blockr)
    result_log['W_FBA'] = sol.objective_value
    result_log['W_FBA_status'] = sol.status
    result_log['W_metabolites'] = len(final_genre.metabolites)
    result_log['W_reactions'] = len(final_genre.reactions)
    result_log['W_genes'] = len(final_genre.genes)
    # save model
    cobra.io.write_sbml_model(final_genre, f'iAF987_{args.prec}_gfW_v0_filterU.xml')

    print(f'__________{args.prec}_____{args.gf}_________')
    print(result_log)
    print('##########################')

### gapfiller original
##########################################################
########### ori pFBA gapfiller from recon_gf.py###########
##########################################################
# draft_genre = model
# universal_obj = str(draft_genre.objective.expression).split()[0].split('*')[-1]

# print('universal_obj:',universal_obj)
elif args.gf == 'ori':
    new_reactions = _find_reactions(draft_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 1, file_type)
    print('@@@ this is new_reactions:',new_reactions)
    result_log['addgf_reactions'] = len(new_reactions)
    inter  = set(rm_rxn_id).intersection(set(new_reactions)) 
    print('intersection:',inter,'total added:',len(new_reactions),'removed reactions:',len(rm_rxn))
    # removed reactions: 21 
    # ['SO4tex', 'FLNDPR2r', 'ACACT7r', '3HAD140', 'AACPS3', 'TAURtex', 'MOGDS', '3OAR120',
    # 'PACCOAL', 'PGSA180', '4HBZtex', 'AACPS7', 'PSD140', '2AGPEAT160', 'BWCOS', '3HAD181', 
    # '3MBtex', 'SULR', '2AGPEAT140', 'OPMEACPD', 'NADH17pp']

    result_log['addgf_reactions_inRM'] = inter
    result_log['addgf_reactions_inRM_num'] = len(inter)
    filled_genre = _gapfill_model(draft_genre, universal, new_reactions, universal_obj, 1)
    print('filled_genre:',len(filled_genre.reactions))

    if file_type != 3:
        print('Identifying new metabolism (Step 2 of 2)...')
        filled_genre = _set_base_inputs(filled_genre, universal)
        media_reactions = _find_reactions(filled_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 2, file_type)
        print('media_reactions:',media_reactions,len(media_reactions))
        final_genre = _gapfill_model(filled_genre, universal, media_reactions, universal_obj, 2)
        print('final_genre:',len(final_genre.reactions))
    else:
        final_genre = _add_annotation(filled_genre, universal_obj)
    alladd = set(new_reactions).union(media_reactions)
    inter = alladd.intersection(set(rm_rxn_id))
    print('alladd:',len(alladd),'intersection:',len(inter))
    result_log['addgf_media_reactions_inRM'] = inter
    result_log['addgf_media_reactions_inter'] = len(inter)    
    sol = final_genre.optimize()
    # blockr = cobra.flux_analysis.find_blocked_reactions(final_genre)
    # result_log['afgf_rm_blocked_reactions'] = len(blockr)
    # print('fba un weighted gapped model:',sol)
    # print('len metabolites:',len(final_genre.metabolites))
    # print('len reactions:',len(final_genre.reactions))
    cobra.io.write_sbml_model(final_genre, f'iAF987_{args.prec}_gf_v0_filterU.xml')
    result_log['gf_FBA'] = sol.objective_value
    result_log['gf_FBA_status'] = sol.status
    result_log['gf_metabolites'] = len(final_genre.metabolites)
    result_log['gf_reactions'] = len(final_genre.reactions)
    result_log['gf_genes'] = len(final_genre.genes)

print(f'__________{args.prec}_____{args.gf}__________')
print(result_log)
print('________________________')

outf = f'result_log_{args.seednum}__test_{args.gf}allunr.txt'
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


