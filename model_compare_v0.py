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

args = parser.parse_args()
print('args:',args.prec,args.seednum)
result_log = {}
# model = cobra.io.read_sbml_model('Ec_core_flux1.xml')
# model.metabolites[:3]
# >>> import cobra
# >>> model = cobra.io.read_sbml_model('/ibex/user/niuk0a/funcarve/cobra/bigg/bigg_universe.xml')
# No objective coefficients in model. Unclear what should be optimized
# >>> print(len(model.reactions))
# 16337
# model = load_model("iYS1720")
# model = load_model("iAB_RBC_283")
# model = load_model("iAF987") #Geobacter metallireducens GS-15
model = cobra.io.read_sbml_model('/ibex/user/niuk0a/funcarve/cobra/iAF987.xml')
solution = model.optimize()
print('#ori model FBA solution:',solution)
result_log['ori_model FBA'] = solution.objective_value
# save 
# cobra.io.write_sbml_model(model, 'iAF987.xml')
# exit()
# print(model.objective)
print(model.medium)
# exit()
## matching genome CP000148.1 
names,seqs = read_fasta("/ibex/user/niuk0a/funcarve/cobra/bigg/aaseq_CP000148.1.txt")
pr2ec , predscore = read_clean_withscore('/ibex/user/niuk0a/funcarve/cobra/aaseq_CP000148.1_maxsep.csv')
with open('biggr2ec.pkl', 'rb') as f:
    biggr2ec = pickle.load(f)
with open('biggec2r.pkl', 'rb') as f:
    biggec2r = pickle.load(f)

# #write to a fasta file
# with open("aaseq_CP000148.1.fasta", "w") as f:
#     for name,seq in zip(names, seqs):
#         f.write(f">{name}\n{seq}\n")
# select unreviewed proteins 
unreviewed_f = '/ibex/user/niuk0a/funcarve/cobra/269799_unreviewed.tsv'
unreviewed = pd.read_csv(unreviewed_f, sep='\t')
print(len(unreviewed))
unreviewed_genes = set(unreviewed['Gene Names'].str.split(' ').sum())
print(len(unreviewed_genes))

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
# seed 1234
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
# result_log['rm reaction'] = [r.id for r in rm_rxn]


# print(f"{r}:{gene}".format(r=r for r in rm_rxn, gene=rxn_genes[r for r in rm_rxn]))
# add ={'NADHPO', 'ASPO6', 'ACOAH', 'HSDA'}
# for r in add:
#     try:
#         print(f"{r}:{genes_rxn[r]}")
#     except:
#         continue
# add1 ={'AAMYL', 'TREY'}
# for r in add1:
#     try:
#         print(f"add1 :{r}:{genes_rxn[r]}")
#     except:
#         continue
# exit()
# for r in rm_rxn:
#     rid = r.id
#     print(f"{rid}:{rxn_genes[rid]}")
# print('rm_rxn in r2ec',len([r for r in rm_rxn if r.id in biggr2ec.keys()]),'total:',len(rm_rxn))
# blockr = cobra.flux_analysis.find_blocked_reactions(model)
# print('bf removeblockr:',blockr)      


# print(model.reactions[:10])
# f = '/ibex/user/niuk0a/funcarve/cobra/SEED2VMH_translation.csv'
# ms2bigg={}
# with open(f, 'r') as file:
#     for line in file:
#         line = line.strip('\n').split(',')
#         ms2bigg[line[0]] = line[1]

# universal = cobra.io.read_sbml_model("/ibex/user/niuk0a/funcarve/cobra/bigg/bigg_universe.xml")
universal = cobra.io.load_json_model("/ibex/user/niuk0a/funcarve/cobra/bigg/universal_model_cobrapy.json")
# filename='/ibex/user/niuk0a/funcarve/cobra/universal_recon.pickle'
# with open(filename, 'rb') as f: reconuniversal = pickle.load(f)

# # filter universal
# print('len(reconuniversal.reactions)',len(reconuniversal.reactions))
# uuu  = len(universal.reactions)
# print('len(universal.reactions)',len(universal.reactions))
# urm = []
# for r in reconuniversal.reactions:
#     msid = r.id
#     print('msid:',msid)
#     try:
#         biggid = ms2bigg[msid]
#         print('biggid:',biggid)
#         ur = universal.reactions.get_by_id(biggid)
#         ur.id = biggid
#     except:
#         urm.append(r)
# print('len(urm)',len(urm))
# universal.remove_reactions(urm)
# universal.repair()
# #save universal
# print('len(urm)',len(urm))
# cobra.io.write_sbml_model(universal, 'filter_bigg_universe.xml')
# print('len(reconuniversal.reactions)',len(reconuniversal.reactions))

# print('len(universal.reactions)','before',uuu,'now',len(universal.reactions))  
    
# exit()


print('len(universal.reactions)',len(universal.reactions))
# universal_scoredict = clean2biggr(universal,predscore)
universal_scoredict = clean2biggr(predscore)
min_old = min(universal_scoredict.values())
max_old = max(universal_scoredict.values())

# c=0
# for r in rm_rxn:
#     if r.id in universal_scoredict:
#         print('r:',r.id,universal_scoredict[r.id],10 * (1 - (universal_scoredict[r.id] - min_old) / (max_old - min_old)))
#         c+=1
# print('c:',c,'total:',len(rm_rxn))
print('rm_rxn in r2ec',len([r for r in rm_rxn if r.id in biggr2ec.keys()]),'total:',len(rm_rxn))

# uni = Model("universal")
# for r in rm_rxn:
#     print('r:',r,r.id,r.gene_reaction_rule)
#     uni.add_reactions([model.reactions.get_by_id(r.id).copy()])
#     break

# print('uni reactions:',len(uni.reactions))

# fvb_result = cobra.flux_analysis.flux_variability_analysis(model,fraction_of_optimum=0.9)
# print('bf remove rxn fvb_result:',fvb_result)
# print('rm_rxn[0]:',rm_rxn[0],type(rm_rxn[0]))
# r = rm_rxn[0]
# uni.add_reactions([r.copy()])
# model.remove_reactions([r])
# print('len model reactions:',len(model.reactions))
# print('len metabolites:',len(model.metabolites))
model.remove_reactions(rm_rxn)
model,rm_metabolits = cobra.manipulation.prune_unused_metabolites(model)
# print('len model reactions:',len(model.reactions))
# print('len metabolites:',len(model.metabolites))
# blockr = cobra.flux_analysis.find_blocked_reactions(model)
# result_log['rm_blocked_reactions'] = len(blockr)
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
    # print('r:',r,'|',r.id,sep='')
    if 'BIOMASS' in r.id:
        print('r|',r.id,sep='')
        model.objective = r
if oriobj == model.objective:
    print('!!no biomass reaction found')
    # exit()
else:
    print('biomass reaction found:',model.objective)

print('model.objective.expression:',model.objective.expression)
model.solver.configuration.tolerances.feasibility = 1e-6
# model.solver.configuration.tolerances.optimality = 1e-6

# model.summary()
# print('model opt summ pass ')
# solution = model.optimize()
# print('solution:',solution,'solution.status:',solution.status)
# status = model.solver.status
# print('model.status:',status)
# blockr = cobra.flux_analysis.find_blocked_reactions(model)
# print('blockr:',blockr)
# fvb_result = cobra.flux_analysis.flux_variability_analysis(model,fraction_of_optimum=0.9)
# print('fvb_result:',fvb_result)
# gfrxn = gapfill(model,uni)
# print('gapfilled reactions:',gfrxn)



# gfrxn = gapfill(model, universal)
# print('gapfilled reactions:',gfrxn)
# pena_gfrxn = gapfill(model, universal, penalties={'universal':3,'GF6PTA': 10},iterations=1)
# print(':gfrxn:',gfrxn)

######################################################
########### pFBA gapfiller from recon_gf.py###########
##############################################
# from recon_gf import _find_reactions, _gapfill_model, _set_base_inputs
# draft_genre = model
# metabolic_tasks = [] # not required
# universal_obj = model.objective
# print('universal_obj:',universal_obj)
# universal_obj = str(draft_genre.objective.expression).split()[0].split('*')[-1]
# print('universal_obj:',universal_obj)   
# min_frac = 0.01
# max_frac = 0.5
# file_type = 3
# print('universal_obj:',universal_obj)
# universal.add_reactions([model.reactions.get_by_id('BIOMASS_Gm_GS15_core_79p20M')])
# print('add biomass to universal')
# new_reactions = _find_reactions(draft_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 1, file_type)
# print(' this is new_reactions:',new_reactions)
# print('________________________')

# filled_genre = _gapfill_model(draft_genre, universal, new_reactions, universal_obj, 1)
# print('filled_genre:',len(filled_genre.reactions))
# print('________________________')
# if file_type != 3:
#     print('Identifying new metabolism (Step 2 of 2)...')
#     filled_genre = _set_base_inputs(filled_genre, universal)
#     media_reactions = _find_reactions(filled_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 2, file_type)
#     print('media_reactions:',media_reactions,len(media_reactions))
#     final_genre = _gapfill_model(filled_genre, universal, media_reactions, universal_obj, 2)
#     print('final_genre:',len(final_genre.reactions))
# exit()
######################################################
########### pFBA gapfiller from recon_gf.py###########
######################################################
######### media ########
# Dictionary for ModelSEED to BiGG ID mapping
# media = 'minimal'
# modelseed_to_bigg = {
#     'cpd00001_e': 'h2o_e',
#     'cpd00035_e': 'o2_e',
#     'cpd00041_e': 'co2_e',
#     'cpd00023_e': 'glc__D_e',
#     'cpd00119_e': 'nh4_e',
#     'cpd00107_e': 'pi_e',
#     'cpd00060_e': 'so4_e',
#     'cpd00161_e': 'k_e',
#     'cpd00069_e': 'fe2_e',
#     'cpd00084_e': 'h_e',
#     'cpd00033_e': 'ac_e',
#     'cpd00322_e': 'cl_e',
#     'cpd00066_e': 'mg2_e',
#     'cpd00054_e': 'na1_e',
#     'cpd00065_e': 'ca2_e',
#     'cpd00156_e': 'cu2_e',
#     'cpd00220_e': 'mn2_e',
#     'cpd00644_e': 'zn2_e',
#     'cpd00393_e': 'cobalt2_e',
#     'cpd00133_e': 'ni2_e',
#     'cpd00263_e': 'mobd_e',
#     'cpd00104_e': 'trp__L_e',
#     'cpd00149_e': 'his__L_e',
#     'cpd00971_e': 'gly_e',
#     'cpd00099_e': 'ala__L_e',
#     'cpd00205_e': 'ser__L_e',
#     'cpd00009_e': 'nad_e',
#     'cpd00063_e': 'asp__L_e',
#     'cpd00254_e': 'glu__L_e',
#     'cpd10515_e': 'phe__L_e',
#     'cpd00030_e': 'arg__L_e',
#     'cpd00242_e': 'leu__L_e',
#     'cpd00226_e': 'ile__L_e',
#     'cpd01242_e': 'val__L_e',
#     'cpd00307_e': 'thr__L_e',
#     'cpd00092_e': 'lys__L_e',
#     'cpd00117_e': 'met__L_e',
#     'cpd00067_e': 'pro__L_e',
#     'cpd00567_e': 'cys__L_e',
#     'cpd00132_e': 'asn__L_e',
#     'cpd00210_e': 'gln__L_e',
#     'cpd00320_e': 'tyr__L_e',
#     'cpd03279_e': 'orn_e',
#     'cpd00246_e': 'tryptamine_e',
#     'cpd00311_e': 'pyr_e',
#     'cpd00367_e': 'glycine_e',
#     'cpd00277_e': 'phenylalanine_e',
#     'cpd00182_e': 'oxaloacetate_e',
#     'cpd00654_e': 'oxalate_e',
#     'cpd00412_e': 'fumarate_e',
#     'cpd00438_e': 'succinate_e',
#     'cpd00274_e': 'malate_e',
#     'cpd00186_e': 'citrate_e',
#     'cpd00637_e': 'lactate_e',
#     'cpd00105_e': 'aspartate_e',
#     'cpd00305_e': 'alpha_KG_e',
#     'cpd00309_e': 'pyruvate_e',
#     'cpd00098_e': 'glutamate_e',
#     'cpd00207_e': 'glutamine_e',
#     'cpd00082_e': 'formate_e',
#     'cpd00129_e': 'succ_e'
# }

# Your media setup code
# if media == 'rich':
#     media = ['cpd00001_e','cpd00035_e','cpd00041_e','cpd00023_e','cpd00119_e','cpd00107_e','cpd00060_e',
#              'cpd00161_e','cpd00069_e','cpd00084_e','cpd00033_e', 'cpd00322_e','cpd00066_e','cpd00054_e',
#              'cpd00065_e','cpd00156_e','cpd00220_e','cpd00644_e','cpd00393_e','cpd00133_e','cpd00263_e',
#              'cpd00104_e','cpd00149_e','cpd00971_e','cpd00099_e','cpd00205_e','cpd00009_e','cpd00063_e',
#              'cpd00254_e','cpd10515_e','cpd00030_e','cpd00242_e','cpd00226_e','cpd01242_e','cpd00307_e',
#              'cpd00092_e','cpd00117_e','cpd00067_e','cpd00567_e','cpd00132_e','cpd00210_e','cpd00320_e',
#              'cpd03279_e','cpd00246_e','cpd00311_e','cpd00367_e','cpd00277_e','cpd00182_e','cpd00654_e',
#              'cpd00412_e','cpd00438_e','cpd00274_e','cpd00186_e','cpd00637_e','cpd00105_e','cpd00305_e',
#              'cpd00309_e','cpd00098_e','cpd00207_e','cpd00082_e','cpd00129_e']
# elif media == 'minimal':
#     media = ['cpd00001_e','cpd00065_e','cpd00060_e','cpd00322_e','cpd00129_e','cpd00156_e','cpd00107_e',
#              'cpd00084_e','cpd00149_e','cpd00099_e','cpd10515_e','cpd00030_e','cpd00254_e','cpd00063_e',
#              'cpd00205_e','cpd00009_e','cpd00971_e','cpd00242_e','cpd00104_e','cpd00644_e','cpd00263_e',
#              'cpd00082_e']
# elif media == 'default':
#     media = ['cpd00035_e','cpd00051_e','cpd00132_e','cpd00041_e','cpd00084_e','cpd00053_e','cpd00023_e',
#              'cpd00033_e','cpd00119_e','cpd00322_e','cpd00107_e','cpd00039_e','cpd00060_e','cpd00066_e',
#              'cpd00129_e','cpd00054_e','cpd00161_e','cpd00065_e','cpd00069_e','cpd00156_e','cpd00027_e',
#              'cpd00149_e','cpd00030_e','cpd00254_e','cpd00971_e','cpd00063_e','cpd10515_e','cpd00205_e',
#              'cpd00099_e']
# else:
#     media = media

# # Replace ModelSEED IDs with BiGG IDs
# media = [modelseed_to_bigg[cpd] if cpd in modelseed_to_bigg else cpd for cpd in media]
# print('media:',media)

# # # Set media condition
# if len(media) != 0:
#     media_condition = set(['EX_' + cpd for cpd in media])
#     universal_reactions = set([x.id for x in universal.reactions])
#     for rxn in universal_reactions:
#         if rxn.startswith('EX_'):
#             universal.reactions.get_by_id(rxn).bounds = (0, 1000.0)
#         if rxn in media_condition:
#             universal.reactions.get_by_id(rxn).bounds = (-1000.0, 10000)
media = model.medium.keys()
######## media ########
######################################################
########### weighted pFBA gapfiller from recon_gf.py###########
##############################################
from recon_gf import _find_reactions, _gapfill_model, _set_base_inputs
draft_genre = model
# save draft model
# cobra.io.write_sbml_model(draft_genre, f'iAF987_{args.seednum}_{args.prec}_rm.xml')
# print('saved draft model')
# exit()
## convert bigg model to Modelseed model



metabolic_tasks = [] # not required
universal_obj = model.objective
universal_obj = str(draft_genre.objective.expression).split()[0].split('*')[-1]
print('universal_obj:',universal_obj)   
min_frac = 0.01
max_frac = 0.5
file_type = 1
universal.add_reactions([model.reactions.get_by_id('BIOMASS_Gm_GS15_core_79p20M')])
new_reactions ,coefficientDict = weighted_find_reactions(draft_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 1, file_type, universal_scoredict,rm_rxn)
print('@@@ this is weighted new_reactions ids:',new_reactions)
result_log['addW_reactions'] = len(new_reactions)
# reaction are removed from model
inter  = set(rm_rxn_id).intersection(set(new_reactions)) 
print('intersection:',inter,'total added:',len(new_reactions),'removed reactions:',len(rm_rxn))
result_log['addW_reactions_inRM'] = inter
result_log['addW_reactions_inRM_num'] = len(inter)

filled_genre = _gapfill_model(draft_genre, universal, new_reactions, universal_obj, 1)
print('filled_genre:',len(filled_genre.reactions))
sol = filled_genre.optimize()
print('fba weighted gapped filled_genre model:',sol.objective_value)

## add gf 
# gfrs = gapfill(filled_genre, universal, penalties=coefficientDict,demand_reactions = False, iterations=1 )
# gfrsid = [r.id for r in gfrs]
# ii = set(gfrsid).intersection(set(rm_rxn_id))
# print('filled_genre1:',len(gfrs))


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
# cobra.io.write_sbml_model(final_genre, f'iAF987_{args.prec}_gfW.xml')
print(f'__________{args.prec}______________')
print(result_log)
print('##########################')

### gapfiller original
draft_genre = model
universal_obj = str(draft_genre.objective.expression).split()[0].split('*')[-1]

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
cobra.io.write_sbml_model(final_genre, f'iAF987_{args.prec}_gf.xml')



print('________________________')
result_log['gf_FBA'] = sol.objective_value
result_log['gf_FBA_status'] = sol.status
result_log['gf_metabolites'] = len(final_genre.metabolites)
result_log['gf_reactions'] = len(final_genre.reactions)
result_log['gf_genes'] = len(final_genre.genes)
outf = f'result_log_{args.seednum}_allunr.txt'
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



# ######################################################
# ########### weighted pFBA gapfiller from recon_gf.py###########
# ######################################################
# # random remove 5% reactions which have gene associated
# total_rxn = len(model.reactions)
# num_rm_rxn = int(total_rxn * 0.05)
# #select 5% reactions randomly
# import random
# rm_rxn = random.sample(model.reactions, num_rm_rxn)
# genes_rxn = {}
# rxn_genes = {}
# genes_rm_set = set()
# for x in rm_rxn:
#     if x.gene_reaction_rule != '':
#         if 'or' in x.gene_reaction_rule or 'and' in x.gene_reaction_rule:
#             if 'and' in x.gene_reaction_rule:
#                 genes = x.gene_reaction_rule.split(' and ')
#             else:
#                 genes = x.gene_reaction_rule.split(' or ')
#             genes_rm_set.update(genes)

#             for gene in genes:
#                 if gene in genes_rxn:
#                     genes_rxn[gene].append(x)
#                 else:
#                     genes_rxn[gene] = [x]
#                 rxn_genes[x] = genes
#         else:
#             genes = x.gene_reaction_rule
#             genes_rxn[genes] = x
#             rxn_genes[x] = genes
#             genes_rm_set.add(genes)
# #remove reactions
# # model.remove_reactions(rm_rxn)
# # print(f"removed {len(rm_rxn)} reactions:{rm_rxn}")
# # print(f"removed {len(genes_rm_set)} genes:{genes_rm_set}")







# # genes
# print("Genes")
# print("-----")
# for x in model.genes:
#     print("%s : %s" % (x.id, x.name))
#     break
# r = model.genes.STM1707.reactions 
# ## find reactions associated with a gene
# print("Reactions")
# print("---------")
# xx =0
# for x in model.reactions:
#     xx +=1
#     print("%s : %s" % (x.id, x.reaction))
    
#     if x.gene_reaction_rule != '':
#         # then the reaction is related to one /more genes
#         if 'or' in x.gene_reaction_rule:
#             genes = x.gene_reaction_rule.split(' or ')
#             proteins = genes_proteins(genes)
            
#     if xx >= 10:
#         break
# # exit()
# rm_RXN = ['GF6PTA','3OAR60']
# r = model.reactions.get_by_id('GF6PTA')
# print('rr;',r)

# for rxn in rm_RXN:
#     r = model.reactions.get_by_id(rxn)
#     model.remove_reactions([r])
#     universal.add_reactions([r.copy()])
# # universal.add_reactions([model.reactions.GF6PTA.copy()])
# # model.remove_reactions([model.reactions.GF6PTA])
# gfrxn = gapfill(model, universal,)
# #cobra.flux_analysis.gapfilling.gapfill(model: cobra.core.Model, 
# #                                       universal: Optional[cobra.core.Model] = None, 
# # lower_bound: float = 0.05, 
# # penalties: Optional[Dict[str, cobra.Reaction]] = None, 
# # demand_reactions: bool = True, 
# # exchange_reactions: bool = False, 
# # iterations: int = 1)[source]¶
# pena_gfrxn = gapfill(model, universal, penalties={'universal':3,'GF6PTA': 10},iterations=2)
# print(':gfrxn:',gfrxn)
# print(':pena_gfrxn:',pena_gfrxn)


# # model = Model('example_model')

# # reaction = Reaction('R_3OAS140')
# # reaction.name = '3 oxoacyl acyl carrier protein synthase n C140 '
# # reaction.subsystem = 'Cell Envelope Biosynthesis'
# # reaction.lower_bound = 0.  # This is the default
# # reaction.upper_bound = 1000.  # This is the default

# # reaction.gene_reaction_rule = '( STM2378 or STM1197 )'
# # reaction.genes
# # print(f'{len(model.reactions)} reactions initially')
# # print(f'{len(model.metabolites)} metabolites initially')
# # print(f'{len(model.genes)} genes initially')

# # model.add_reactions([reaction])

# # # The objects have been added to the model
# # print(f'{len(model.reactions)} reactions')
# # print(f'{len(model.metabolites)} metabolites')
# # print(f'{len(model.genes)} genes')

# # print("Reactions")
# # print("---------")
# # for x in model.reactions:
# #     print("%s : %s" % (x.id, x.reaction))

# # print("")
# # print("Metabolites")
# # print("-----------")
# # for x in model.metabolites:
# #     print('%9s : %s' % (x.id, x.formula))

# # print("")
# # print("Genes")
# # print("-----")
# # for x in model.genes:
# #     associated_ids = (i.id for i in x.reactions)
# #     print("%s is associated with reactions: %s" %
# #           (x.id, "{" + ", ".join(associated_ids) + "}"))

# # model.objective = 'R_3OAS140'
# # print(model.objective.expression)
# # print(model.objective.direction)

# # import tempfile
# # from pprint import pprint
# # from cobra.io import write_sbml_model, validate_sbml_model
# # with tempfile.NamedTemporaryFile(suffix='.xml') as f_sbml:
# #     write_sbml_model(model, filename=f_sbml.name)
# #     report = validate_sbml_model(filename=f_sbml.name)

# # pprint(report)

# # solution = model.optimize()
# # print(solution)

# # print(model.objective.expression)
# # print(model.objective.direction)

# # import tempfile
# # from pprint import pprint
# # from cobra.io import write_sbml_model, validate_sbml_model
# # with tempfile.NamedTemporaryFile(suffix='.xml') as f_sbml:
# #     write_sbml_model(model, filename=f_sbml.name)
# #     report = validate_sbml_model(filename=f_sbml.name)

# # pprint(report)

# # solution = model.optimize()
# # print(solution)

