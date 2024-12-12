    
from copy import deepcopy
import numpy as np
import pickle
from sys import stdout
import warnings
import pandas as pd
import symengine
import cobra
# from CLEAN.utils import *
# from CLEAN.infer import infer_maxsep
from cobra.manipulation.delete import *
import os

data_root = '../data/'

with open(f'{data_root}biggr2ec.pkl', 'rb') as f:
    biggr2ec = pickle.load(f)
with open(f'{data_root}biggec2r.pkl', 'rb') as f:
    biggec2r = pickle.load(f)
print("bigg-ec-reaction dictionary loaded successfully.")
with open(f'{data_root}seedr2ec.pkl', 'rb') as f:
    seedr2ec = pickle.load(f)
with open(f'{data_root}seedec2r.pkl', 'rb') as f:
    seedec2r = pickle.load(f)
print("seed-ec-reaction dictionary loaded successfully.")

def run_clean(fasta_file):
    # print('Running CLEAN...',flush=True)
    print('running CLEAN AUTO IS NOT DONE YET',flush=True)
    pklfileofcleanresult= ''
    # fastaname = fasta_file.split('/')[-1].split('.')[0]
    # prepare_infer_fasta(fastaname)
    # infer_maxsep('split100', fastaname, report_metrics=False, pretrained=True, gmm = './data/pretrained/gmm_ensumble.pkl')
    # os.remove("data/"+ fastaname +'.csv')
    # cmd_line = 'python /ibex/user/niuk0a/anaconda3/envs/recon/lib/python3.9/site-packages/reconstructor/CLEAN.py --input ' + fasta_file + ' --output ' + fasta_file.rstrip('.faa') + '_maxsep.csv'
    # os.system(cmd_line)
    # fasta_file.rstrip('.faa') + '_maxsep.csv'
    return pklfileofcleanresult

def read_clean_withscore(input_file,threshold=0.5):
    print('threrhold-->',threshold)
    pr2ec = {}
    allec2pr = {}
    ec2pr = {}
    allpr2ec = {}
    with open(input_file, 'r') as inFile:
        for line in inFile:
            line = line.strip('\n')
            line = line.split(',')
            pr = line[0]
            # items = line[-1].split(',')
            items = line[1:]

            for item in items:
                if item.startswith('EC:'):
                    ec,dis = item.split('/')
                    ecid = ec.split(':')[-1]
                    # ecid = ec
                    dis = float(dis)
                    try:
                        allpr2ec[pr].update({ecid:dis})
                    except:
                        allpr2ec[pr] = {ecid:dis}
                    try:
                        allec2pr[ecid].update({pr:dis})
                    except:
                        allec2pr[ecid] = {pr:dis}

                    if dis <= threshold:
                        try:
                            ec2pr[ecid].update({pr:dis})
                        except:
                            ec2pr[ecid] = {pr:dis}
                        try:
                            pr2ec[pr].append(ecid)
                        except:
                            pr2ec[pr] = [ecid]
    print('pr2ec-protein number->',len(list(pr2ec.keys())))   
    return pr2ec,ec2pr,allpr2ec,allec2pr

def read_cleandf_withscore(input_file,threshold=8):
    df = pd.read_pickle(input_file)

    allpr2ec = df.to_dict()
    allec2pr = df.T.to_dict()
    
    df = df[df <= threshold]
    # df_cleaned = df.dropna(how='all', axis=0).dropna(how='all', axis=1)

    pr2ec = {col: {row: value for row, value in df[col].dropna().items()} 
                                  for col in df.columns
                                  if not df[col].dropna().empty
                                  }
    ec2pr = {row: {col: value for col, value in df.loc[row].dropna().items()} 
                        for row in df.index 
                        if not df.loc[row].dropna().empty}
    
    return pr2ec,ec2pr,allpr2ec,allec2pr

def maximum_separation(dist_lst, first_grad, use_max_grad):
    #cite from CLEAN
    opt = 0 if first_grad else -1
    gamma = np.append(dist_lst[1:], np.repeat(dist_lst[-1], 10))
    sep_lst = np.abs(dist_lst - np.mean(gamma))
    sep_grad = np.abs(sep_lst[:-1]-sep_lst[1:])
    if use_max_grad:
        max_sep_i = np.argmax(sep_grad)
    else:
        large_grads = np.where(sep_grad > np.mean(sep_grad))
        max_sep_i = large_grads[-1][opt]
    if max_sep_i >= 5:
        max_sep_i = 0
    return max_sep_i

def clean2seedr(allpr2ec,threshold):
    universal_scoredict={}
    rxns = {}
    for pr, ec_score in allpr2ec.items():
        smallest_10_dist_df = pd.Series(ec_score).nsmallest(20)
        dist_lst = list(smallest_10_dist_df)
        max_sep_i = maximum_separation(dist_lst, True, True)
        max_sep_ec = [smallest_10_dist_df.index[i] for i in range(max_sep_i+1)]
        for ec , score in ec_score.items():
            ec = ec.split(':')[-1]
            if ec in seedec2r.keys():
                rids = seedec2r[ec]
                for rid in rids:
                    rid = rid + '_c'
                    if rid in universal_scoredict.keys():
                        ori = universal_scoredict[rid]
                        universal_scoredict[rid] = min(score,ori)
                        if universal_scoredict[rid] <= threshold:
                            if ori ==  universal_scoredict[rid]:
                                continue
                            elif ec in max_sep_ec:
                                try:
                                    rxns[rid].append(pr) 
                                except KeyError:
                                    rxns[rid] = [pr]
                    else:
                        universal_scoredict[rid] = score
                        if score <= threshold and ec in max_sep_ec:
                            try:
                                rxns[rid].append(pr) 
                            except KeyError:
                                rxns[rid] = [pr]
    print('rxns:',len(rxns))
    return universal_scoredict,rxns

def update_predscore(newreactions,allec2pr,newallpr2ec,updateprs,reward,threshold):
    count = 0
    for r in newreactions:
        r = r.split('_')[0]
        try:
            ecid = seedr2ec[r]
        except KeyError:
            count+=1
            continue
        for ec in ecid:
            try:
                prd = allec2pr[ec]
            except KeyError:
                count+=1
                continue
            try:
                pr = min(prd, key=prd.get) # take the pr with MIN val for PREDscore
            except ValueError:
                count+=1
                continue
            try:
                s = max(float(newallpr2ec[pr][ec]) - reward, 0)
                if s <= threshold:
                    updateprs.append(pr)
                # print(pr,ec,newallpr2ec[pr][ec],'->', s,flush=True)
                newallpr2ec[pr].update({ec: s})
                allec2pr[ec].update({pr: s})
            except KeyError:
                count+=1
    # print(f'during update: {count} missing')    
    return newallpr2ec,updateprs,allec2pr

def update_predscore_flux(rxn_flux,allec2pr,newpredscore,updateprs,reward,threshold,fluxflage):
    count = 0
    turnovercount = 0
    upp=set()

    if fluxflage == 1:
        print('gradient reward by flux...')
        print('update reactions based on flux ...',flush=True)
    else:
        print('reward is constant val not effect by flux...',flush=True)
    for r in rxn_flux.keys():
        flux = rxn_flux[r]
        r = r.split('_')[0]
        if fluxflage == 1:
            if flux > 0:
                reward = reward*10
            elif flux > 1e-5:
                reward = reward*5
            elif flux > 1e-6:
                reward = reward
        try:
            ecid = seedr2ec[r]
        except KeyError:
            count+=1
            continue
        for ec in ecid:
            try:
                prd = allec2pr[ec]

            except KeyError:
                count+=1
                continue
            try:
                pr = min(prd, key=prd.get)
            except ValueError:
                count+=1
                continue
            try:
                oriscore = float(newpredscore[pr][ec])
                if (oriscore - reward) < threshold:
                    updateprs.append(pr)
                    turnovercount+=1
                s = max(float(newpredscore[pr][ec]) - reward, 0)
                newpredscore[pr].update({ec: s})
                upp.add(pr)
                allec2pr[ec].update({pr: s})
            except KeyError:
                count+=1
    print(f'during update: {count} missing')    
    print(f'In total {turnovercount} proteins are turn overed,{len(upp)} proteins are updated',flush=True)
    return newpredscore,updateprs,allec2pr

def update_predscore_block(blockrxn,allec2pr,newpredscore,updateprs,reward,threshold):
    print('updata block reactions...',flush=True)
    count = 0
    turnovercount = 0
    for r in blockrxn:
        r = r.split('_')[0]
        try:
            ecid = seedr2ec[r]
        except KeyError:
            count+=1
            continue
        for ec in ecid:
            try:
                prd = allec2pr[ec]
            except KeyError:
                count+=1
                continue
            try:
                pr = min(prd, key=prd.get)
            except ValueError:
                count+=1
                continue
            try:
                if float(newpredscore[pr][ec]) + reward > threshold:
                    updateprs.append(pr)
                    turnovercount+=1
                s = float(newpredscore[pr][ec]) + reward
                # print(pr,ec,newpredscore[pr][ec],'->', s,flush=True)
                newpredscore[pr].update({ec: s})
                allec2pr[ec].update({pr: s})

            except KeyError:
                count+=1
    print(f'during update: {count} missing')    
    print(f'In total {turnovercount} proteins are turn overed, removed from the draft model ',flush=True)
    return newpredscore,updateprs,allec2pr
    
def get_org_rxns(gene_modelseed, organism):
    ''' Get genes for organism from reference genome '''
    rxn_db = {}
    org_genes = []
    for gene in gene_modelseed.keys():
        current = gene.split(':')[0]
        if current == organism:
            org_genes.append(gene)

    return set(org_genes)

def create_model(rxn_db, universal, input_id):
    ''' Create draft GENRE and integrate GPRs '''
    new_model = cobra.Model('new_model')
    tmpuniversal = deepcopy(universal)
    for x in rxn_db.keys():
        try:
            rxn = tmpuniversal.reactions.get_by_id(x)
            rxn.gene_reaction_rule = ' or '.join(rxn_db[x])
            new_model.add_reactions([rxn])
        except KeyError:
            continue
    if input_id != 'default':
        new_model.id = input_id
    return new_model

# Add gene names
def add_names(model, gene_db):
    ''' Add gene names '''
    for gene in model.genes:
        try:
            gene.name = gene_db[gene.id].title()
        except KeyError:
            continue
    return model

# pFBA gapfiller
def find_reactions(model, reaction_bag, tasks, obj, fraction, max_fraction, step, file_type):
    ''' pFBA gapfiller that modifies universal reaction bag, removes overlapping reacitons from universal reaction bag
    and resets objective if needed, adds model reaction to universal bag, sets lower bound for metabolic tasks, 
    sets minimum lower bound for previous objective, assemble forward and reverse components of all reactions,
    create objective, based on pFBA, run FBA and identify reactions from universal that are now active'''
    stdout.write('\r[                                         ]')
    stdout.flush()

    # Modify universal reaction bag
    new_rxn_ids = set() #make empty set we will add new reaction ids to
    with reaction_bag as universal: #set the reaction bag as the universal reaction databse

        # Remove overlapping reactions from universal bag, and reset objective if needed
        warnings.filterwarnings('ignore')
        orig_rxn_ids = set()    #original reaction ids start as an empty set
        remove_rxns = []    #reactions to remove is empty vector
        for rxn in model.reactions: #for a reaction in the draft model reactions
            if rxn.id == obj and file_type != 3: #if a reaction is part of the objective function 
                continue

            orig_rxn_ids |= set([rxn.id])
            try:
                test = universal.reactions.get_by_id(rxn.id)
                remove_rxns.append(rxn.id)
            except:
                continue

        # Add model reactions to universal bag
        universal.remove_reactions(list(set(remove_rxns)))
        add_rxns = []
        for x in model.reactions:
            if x.id != obj or file_type == 3:
                add_rxns.append(deepcopy(x))
        universal.add_reactions(add_rxns)

        # Set lower bounds for metaboloic tasks
        if len(tasks) != 0:
            for rxn in tasks:
                try:
                    universal.reactions.get_by_id(rxn).lower_bound = fraction
                except:
                    continue

        stdout.write('\r[---------------                          ]')
        stdout.flush()

        # Set minimum lower bound for previous objective
        universal.objective = universal.reactions.get_by_id(obj) 
        prev_obj_val = universal.slim_optimize()
        print('\nuniversal.objective:', universal.objective)
        print('universal.slim_optimize:', prev_obj_val,flush=True)
        if step == 1:
            print('1 current constraint',universal.reactions.get_by_id(obj))
            prev_obj_constraint = universal.problem.Constraint(universal.reactions.get_by_id(obj).flux_expression, 
        	   lb=prev_obj_val*fraction, ub=prev_obj_val*max_fraction)
        elif step == 2:
            print('2 current constraint',universal.reactions.get_by_id(obj))
            
            prev_obj_constraint = universal.problem.Constraint(universal.reactions.get_by_id(obj).flux_expression, 
               lb=prev_obj_val*max_fraction, ub=prev_obj_val)
        universal.solver.add(prev_obj_constraint)
        universal.solver.update()

        # Assemble forward and reverse components of all reactions
        coefficientDict = {}
        pfba_expr = symengine.RealDouble(0)
        for rxn in universal.reactions:
            if rxn.id in orig_rxn_ids:
                coefficientDict[rxn.forward_variable] = 0.0
                coefficientDict[rxn.reverse_variable] = 0.0
            else:
                coefficientDict[rxn.forward_variable] = 1.0
                coefficientDict[rxn.reverse_variable] = 1.0

        stdout.write('\r[--------------------------               ]')
        stdout.flush()

        # Create objective, based on pFBA
        universal.objective = 0
        universal.solver.update()
        universal.objective = universal.problem.Objective(symengine.RealDouble(0), direction='min', sloppy=True)
        universal.objective.set_linear_coefficients(coefficientDict)
        
        stdout.write('\r[----------------------------------       ]')
        stdout.flush()

        # Run FBA and identify reactions from universal that are now active
        solution = universal.optimize()
        
    # new_rxn_ids = set([rxn.id for rxn in reaction_bag.reactions if abs(solution.fluxes[rxn.id]) > 1e-6]).difference(orig_rxn_ids)
    print('step {}-'.format(step),solution.fluxes[rxn.id])
    print('step {}-'.format(step),len(orig_rxn_ids),list(orig_rxn_ids)[:10])
    coutn=0
    for rxn in reaction_bag.reactions:
        a = solution.fluxes[rxn.id]
        if abs(a) > 1e-13:
            print(rxn.id,'---->',a,end='\n')
        else:
            if a == 0.0:
                coutn+=1
            
    print('step {}-'.format(step),len(reaction_bag.reactions))
    print('step {}-'.format(step),coutn)
    threshold1 = 1e-10
    new_rxn_ids = set([rxn.id for rxn in reaction_bag.reactions if abs(solution.fluxes[rxn.id]) > threshold1]).difference(orig_rxn_ids)
    print('threshold:',threshold1)
    stdout.write('\r[-----------------------------------------]\n')
    warnings.filterwarnings('default')
    print('step {} --new_rxn_ids:==>'.format(step),len(new_rxn_ids),'\n',new_rxn_ids)
    return(new_rxn_ids)     

# pFBA gapfiller
def weighted_find_reactions(model, reaction_bag, tasks, obj, fraction, max_fraction, step, file_type,weight_dict,upper=15,lower=5,maxweight=100,minweight=0.0):
    ''' pFBA gapfiller that modifies universal reaction bag, removes overlapping reacitons from universal reaction bag
    and resets objective if needed, adds model reaction to universal bag, sets lower bound for metabolic tasks, 
    sets minimum lower bound for previous objective, assemble forward and reverse components of all reactions,
    create objective, based on pFBA, run FBA and identify reactions from universal that are now active'''

    # Modify universal reaction bag
    new_rxn_ids = set() #make empty set we will add new reaction ids to
    with reaction_bag as universal: #set the reaction bag as the universal reaction databse

        # Remove overlapping reactions from universal bag, and reset objective if needed
        print('#Remove overlapping reactions from universal bag, and reset objective if needed')
        warnings.filterwarnings('ignore')
        orig_rxn_ids = set()    #original reaction ids start as an empty set
        remove_rxns = []    #reactions to remove is empty vector
        for rxn in model.reactions: #for a reaction in the draft model reactions
            if rxn.id == obj and file_type != 3: #if a reaction is part of the objective function 
                continue
            orig_rxn_ids |= set([rxn.id])
            try:
                test = universal.reactions.get_by_id(rxn.id)
                remove_rxns.append(rxn.id)
            except:
                continue
        # Add model reactions to universal bag
        print('#Add model reactions to universal bag')
        universal.remove_reactions(list(set(remove_rxns)))
        add_rxns = []
        for x in model.reactions:
            if x.id != obj or file_type == 3:
                add_rxns.append(x.copy())
        universal.add_reactions(add_rxns)
        # Set lower bounds for metaboloic tasks
        print('#Set lower bounds for metaboloic tasks')
        if len(tasks) != 0:
            for rxn in tasks:
                try:
                    universal.reactions.get_by_id(rxn).lower_bound = fraction
                except:
                    continue

        # Set minimum lower bound for previous objective
        print('#Set minimum lower bound for previous objective')
        universal.objective = universal.reactions.get_by_id(obj) 
        prev_obj_val = universal.slim_optimize()
        print('#Prev_obj_val:',prev_obj_val)
        # if prev_obj_val < 0:
        #     prev_obj_val = abs(prev_obj_val)
        
        if step == 1:
            prev_obj_constraint = universal.problem.Constraint(universal.reactions.get_by_id(obj).flux_expression, 
        	   lb=prev_obj_val*fraction, ub=prev_obj_val*max_fraction)
        elif step == 2:
            prev_obj_constraint = universal.problem.Constraint(universal.reactions.get_by_id(obj).flux_expression, 
               lb=prev_obj_val*max_fraction, ub=prev_obj_val)
        universal.solver.add(prev_obj_constraint)
        universal.solver.update()
        # Assemble forward and reverse components of all reactions
        print('#Assemble forward and reverse components of all reactions')
        coefficientDict = {}
# upper=15,lower=5,maxweight=100,minweight=0.0 
        for rxn in universal.reactions:
            if rxn.id in orig_rxn_ids:
                coefficientDict[rxn.forward_variable] = minweight
                coefficientDict[rxn.reverse_variable] = minweight
            else:
                if rxn.id in weight_dict:
                    if weight_dict[rxn.id] < lower:
                        coefficientDict[rxn.forward_variable] = minweight
                        coefficientDict[rxn.reverse_variable] = minweight
                    elif weight_dict[rxn.id] < upper:
                        coefficientDict[rxn.forward_variable] = (weight_dict[rxn.id]-lower)*(maxweight-minweight)/(upper-lower)+minweight
                        coefficientDict[rxn.reverse_variable] = (weight_dict[rxn.id]-lower)*(maxweight-minweight)/(upper-lower)+minweight
                else:
                    coefficientDict[rxn.forward_variable] = maxweight
                    coefficientDict[rxn.reverse_variable] = maxweight

        # Create objective, based on pFBA
        print('#Create objective, based on pFBA')
        universal.objective = 0
        universal.solver.update()
        universal.objective = universal.problem.Objective(symengine.RealDouble(0), direction='min', sloppy=True)
        universal.objective.set_linear_coefficients(coefficientDict)
        # Run FBA and identify reactions from universal that are now active
        print('#Run FBA and identify reactions from universal that are now active')
        solution = universal.optimize()
        print('#Run FBA solution:',solution)

    fluxthreshold = 1e-6
    print('Fluxthreshold:',fluxthreshold)
    new_rxn_flux = {}
    new_rxn_ids = set()
    for rxn in reaction_bag.reactions:
        if abs(solution.fluxes[rxn.id]) > fluxthreshold:
            if rxn.id not in orig_rxn_ids:
                new_rxn_ids.add(rxn.id)
                new_rxn_flux[rxn.id] = solution.fluxes[rxn.id]
    warnings.filterwarnings('default')
    print('Return with flux:',len(new_rxn_flux.keys()))
    return new_rxn_ids ,coefficientDict,new_rxn_flux


# Add new reactions to model
def gapfill_model(model, universal, new_rxn_ids, obj, step):
    '''Adds new reactions to model by getting reactions and metabolites to be added to the model, creates gapfilled model, 
    and identifies extracellular metabolites that still need exchanges '''
    # Get reactions and metabolites to be added to the model
    new_rxns = []
    if step == 1: new_rxns.append(deepcopy(universal.reactions.get_by_id(obj)))
    for rxn in new_rxn_ids: 
        if rxn != obj:
            new_rxns.append(deepcopy(universal.reactions.get_by_id(rxn)))
    
    # Create gapfilled model 
    model.add_reactions(new_rxns)
    model.objective = model.problem.Objective(model.reactions.get_by_id(obj).flux_expression, direction='max')

    # Identify extracellular metabolites still need exchanges
    for cpd in model.metabolites:
        if cpd.compartment != 'extracellular':
            continue
        else:
            try:
                test = model.reactions.get_by_id('EX_' + cpd.id)
            except KeyError:
                exch_id = 'EX_' + cpd.id
                model.add_boundary(cpd, type='exchange', reaction_id=exch_id, lb=-1000.0, ub=1000.0)
                model.reactions.get_by_id(exch_id).name = cpd.name + ' exchange'

    return model


# Set uptake of specific metabolites in complete medium gap-filling
def set_base_inputs(model, universal):
    tasks = ['EX_cpd00035_e','EX_cpd00051_e','EX_cpd00132_e','EX_cpd00041_e','EX_cpd00084_e','EX_cpd00053_e','EX_cpd00023_e',
    'EX_cpd00033_e','EX_cpd00119_e','EX_cpd00322_e','EX_cpd00107_e','EX_cpd00039_e','EX_cpd00060_e','EX_cpd00066_e','EX_cpd00129_e',
    'EX_cpd00054_e','EX_cpd00161_e','EX_cpd00065_e','EX_cpd00069_e','EX_cpd00156_e','EX_cpd00027_e','EX_cpd00149_e','EX_cpd00030_e',
    'EX_cpd00254_e','EX_cpd00971_e','EX_cpd00063_e','EX_cpd10515_e','EX_cpd00205_e','EX_cpd00099_e']
    new_rxns = []
    for exch in tasks: 
        try:
            test = model.reactions.get_by_id(exch)
        except:
            new_rxns.append(deepcopy(universal.reactions.get_by_id(exch)))
    model.add_reactions(new_rxns)
    for exch in tasks: model.reactions.get_by_id(exch).bounds = (-1000., -0.01)
    return model

def add_annotation(model, gram, obj='built'):
    ''' Add gene, metabolite, reaction ,biomass reaction annotations '''
    # Genes
    for gene in model.genes:
        gene._annotation = {}
        gene.annotation['sbo'] = 'SBO:0000243'
        gene.annotation['kegg.genes'] = gene.id
    
    # Metabolites
    for cpd in model.metabolites: 
        cpd._annotation = {}
        cpd.annotation['sbo'] = 'SBO:0000247'
        if 'cpd' in cpd.id: cpd.annotation['seed.compound'] = cpd.id.split('_')[0]

    # Reactions
    for rxn in model.reactions:
        rxn._annotation = {}
        if 'rxn' in rxn.id: rxn.annotation['seed.reaction'] = rxn.id.split('_')[0]
        compartments = set([x.compartment for x in list(rxn.metabolites)])
        if len(list(rxn.metabolites)) == 1:
            rxn.annotation['sbo'] = 'SBO:0000627' # exchange
        elif len(compartments) > 1:
            rxn.annotation['sbo'] = 'SBO:0000185' # transport
        else:
            rxn.annotation['sbo'] = 'SBO:0000176' # metabolic

    # Biomass reactions
    if obj == 'built':
        try:
            model.reactions.EX_biomass.annotation['sbo'] = 'SBO:0000632'
        except:
            pass
        if gram == 'none':
            biomass_ids = ['dna_rxn','rna_rxn','protein_rxn','teichoicacid_rxn','lipid_rxn','cofactor_rxn','rxn10088_c','biomass_rxn']
        else:
            biomass_ids = ['dna_rxn','rna_rxn','protein_rxn','teichoicacid_rxn','peptidoglycan_rxn','lipid_rxn','cofactor_rxn','GmPos_cellwall','rxn10088_c','GmNeg_cellwall','biomass_rxn_gp','biomass_rxn_gn']
        for x in biomass_ids:
            try:
                model.reactions.get_by_id(x).annotation['sbo'] = 'SBO:0000629'
            except:
                continue
    else:
        model.reactions.get_by_id(obj).annotation['sbo'] = 'SBO:0000629'

    return model
    

# Run some basic checks on new models
def checkModel(pre_reactions, pre_metabolites, post_model):
    print('\n\tChecking new model...',flush=True)  
    ''' Run basic checks on new models (checking for objective flux'''

	# Check for objective flux
    new_genes = len(post_model.genes)
    new_rxn_ids = set([x.id for x in post_model.reactions]).difference(pre_reactions)
    new_cpd_ids = set([x.id for x in post_model.metabolites]).difference(pre_metabolites)
    test_flux = round(post_model.slim_optimize(), 3)

	
	# Report to user
    print('\tDraft reconstruction had', str(new_genes), 'genes,', str(len(pre_reactions)), 'reactions, and', str(len(pre_metabolites)), 'metabolites')
    print('\tGapfilled', str(len(new_rxn_ids)), 'reactions and', str(len(new_cpd_ids)), 'metabolites\n')
    print('\tFinal reconstruction has', str(len(post_model.reactions)), 'reactions and', str(len(post_model.metabolites)), 'metabolites')
    print('\tFinal objective flux is', str(round(test_flux, 3)))

    return test_flux

def print_model_info(model):
    print('*'*10)
    print('modle name:',model.name)
    print('number of reactions:',len(model.reactions))
    print('number of metabolites:',len(model.metabolites))
    print('number of genes:',len(model.genes))
    print('objective:',str(model.objective)[:100])
    print('objective expression:',str(model.objective.expression)[:100])
    # print('blocked reactions:',len(cobra.flux_analysis.variability.find_blocked_reactions(model)))
    print('status:',model.solver.status)
    print('*'*10,flush=True)

def differ(predscore,newpredscore):
    differdict={}
    for pr in predscore.keys():
        for ec in predscore[pr].keys():
            try:
                if predscore[pr][ec] != newpredscore[pr][ec]:
                    differdict[pr] = {ec:predscore[pr][ec]}
            except KeyError:
                continue
    return differdict

def infer_confidence_gmm(distance, gmm_lst):
    confidence = []
    for j in range(len(gmm_lst)):
        main_GMM = gmm_lst[j]
        a, b = main_GMM.means_
        true_model_index = 0 if a[0] < b[0] else 1
        certainty = main_GMM.predict_proba([[distance]])[0][true_model_index]
        confidence.append(certainty)
    return np.mean(confidence)

def convert_distance_to_confidence(df,gmmf,first_grad=True, use_max_grad=False):
    '''Convert distance to confidence'''
    for col in df.columns:
        ec=[]
        smallest_10_dist_df = df[col].nsmallest(20)
        dist_lst = list(smallest_10_dist_df)
        max_sep_i = maximum_separation(dist_lst, first_grad, use_max_grad)
        for i in range(max_sep_i+1):
            EC_i = smallest_10_dist_df.index[i]
            dist_i = smallest_10_dist_df[i]
            if gmmf != None:
                gmm_lst = pickle.load(open(gmmf, 'rb'))
                dist_i = infer_confidence_gmm(dist_i, gmm_lst)
            dist_str = "{:.4f}".format(dist_i)