import warnings
from copy import deepcopy
# from cobra.core import Reaction, Metabolite
import cobra
import symengine
from cobra.util import solver
from cobra.flux_analysis import gapfill
from cobra.manipulation.delete import *
import pickle

# Dictionary for ModelSEED to BiGG ID mapping


# # # Set media condition
# # if len(media) != 0:
# #     media_condition = set(['EX_' + cpd for cpd in media])
# #     universal_reactions = set([x.id for x in universal.reactions])
# #     for rxn in universal_reactions:
# #         if rxn.startswith('EX_'):
# #             universal.reactions.get_by_id(rxn).bounds = (0, 1000.0)
# #         if rxn in media_condition:
# #             universal.reactions.get_by_id(rxn).bounds = (-1000.0, 10000)

# pFBA gapfiller
def _find_reactions(model, reaction_bag, tasks, obj, fraction, max_fraction, step, file_type):
    ''' pFBA gapfiller that modifies universal reaction bag, removes overlapping reacitons from universal reaction bag
    and resets objective if needed, adds model reaction to universal bag, sets lower bound for metabolic tasks, 
    sets minimum lower bound for previous objective, assemble forward and reverse components of all reactions,
    create objective, based on pFBA, run FBA and identify reactions from universal that are now active'''
    # stdout.write('\r[                                         ]')
    # stdout.flush()

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
        print('#add model reactions to universal bag')
        universal.remove_reactions(list(set(remove_rxns)))
        add_rxns = []
        for x in model.reactions:
            if x.id != obj or file_type == 3:
                add_rxns.append(deepcopy(x))
        universal.add_reactions(add_rxns)

        # Set lower bounds for metaboloic tasks
        print('#set lower bounds for metaboloic tasks')
        if len(tasks) != 0:
            for rxn in tasks:
                try:
                    universal.reactions.get_by_id(rxn).lower_bound = fraction
                except:
                    continue

        # stdout.write('\r[---------------                          ]')
        # stdout.flush()

        # Set minimum lower bound for previous objective
        print('#set minimum lower bound for previous objective')
        universal.objective = universal.reactions.get_by_id(obj) 
        prev_obj_val = universal.slim_optimize()
        print('#prev_obj_val:',prev_obj_val)
        if step == 1:
            prev_obj_constraint = universal.problem.Constraint(universal.reactions.get_by_id(obj).flux_expression, 
        	   lb=prev_obj_val*fraction, ub=prev_obj_val*max_fraction)
        elif step == 2:
            prev_obj_constraint = universal.problem.Constraint(universal.reactions.get_by_id(obj).flux_expression, 
               lb=prev_obj_val*max_fraction, ub=prev_obj_val)
        universal.solver.add(prev_obj_constraint)
        universal.solver.update()

        # Assemble forward and reverse components of all reactions
        print('#assemble forward and reverse components of all reactions')
        coefficientDict = {}
        pfba_expr = symengine.RealDouble(0)
        for rxn in universal.reactions:
            if rxn.id in orig_rxn_ids:
                coefficientDict[rxn.forward_variable] = 0.0
                coefficientDict[rxn.reverse_variable] = 0.0
            else:
                coefficientDict[rxn.forward_variable] = 1.0
                coefficientDict[rxn.reverse_variable] = 1.0

        # stdout.write('\r[--------------------------               ]')
        # stdout.flush()

        # Create objective, based on pFBA
        print('#create objective, based on pFBA')
        universal.objective = 0
        universal.solver.update()
        universal.objective = universal.problem.Objective(symengine.RealDouble(0), direction='min', sloppy=True)
        universal.objective.set_linear_coefficients(coefficientDict)
        
        # stdout.write('\r[----------------------------------       ]')
        # stdout.flush()

        # Run FBA and identify reactions from universal that are now active
        print('#run FBA and identify reactions from universal that are now active')
        solution = universal.optimize()
        print('# run FBA solution:',solution)
    new_rxn_ids = set([rxn.id for rxn in reaction_bag.reactions if abs(solution.fluxes[rxn.id]) > 1e-6]).difference(orig_rxn_ids)
    # new_rxn_ids = set([rxn.id for rxn in reaction_bag.reactions if abs(solution.fluxes[rxn.id]) > 1e-10]).difference(orig_rxn_ids)
    # stdout.write('\r[-----------------------------------------]\n')
    warnings.filterwarnings('default')

    return(new_rxn_ids)    

# pFBA gapfiller
def weighted_find_reactions(model, reaction_bag, tasks, obj, fraction, max_fraction, step, file_type,weight_dict,rm_rxn):
    ''' pFBA gapfiller that modifies universal reaction bag, removes overlapping reacitons from universal reaction bag
    and resets objective if needed, adds model reaction to universal bag, sets lower bound for metabolic tasks, 
    sets minimum lower bound for previous objective, assemble forward and reverse components of all reactions,
    create objective, based on pFBA, run FBA and identify reactions from universal that are now active'''
    # stdout.write('\r[                                         ]')
    # stdout.flush()

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
        print('#add model reactions to universal bag')
        universal.remove_reactions(list(set(remove_rxns)))
        add_rxns = []
        for x in model.reactions:
            if x.id != obj or file_type == 3:

                add_rxns.append(x.copy())
        universal.add_reactions(add_rxns)
      


        # Set lower bounds for metaboloic tasks
        print('#set lower bounds for metaboloic tasks')
        if len(tasks) != 0:
            for rxn in tasks:
                try:
                    universal.reactions.get_by_id(rxn).lower_bound = fraction
                except:
                    continue

        # stdout.write('\r[---------------                          ]')
        # stdout.flush()

        # Set minimum lower bound for previous objective
        print('#set minimum lower bound for previous objective')
        universal.objective = universal.reactions.get_by_id(obj) 
        prev_obj_val = universal.slim_optimize()
        print('#prev_obj_val:',prev_obj_val)
        if step == 1:
            prev_obj_constraint = universal.problem.Constraint(universal.reactions.get_by_id(obj).flux_expression, 
        	   lb=prev_obj_val*fraction, ub=prev_obj_val*max_fraction)
        elif step == 2:
            prev_obj_constraint = universal.problem.Constraint(universal.reactions.get_by_id(obj).flux_expression, 
               lb=prev_obj_val*max_fraction, ub=prev_obj_val)
        universal.solver.add(prev_obj_constraint)
        universal.solver.update()

        # Assemble forward and reverse components of all reactions
        print('#assemble forward and reverse components of all reactions')
        coefficientDict = {}
        # with open('biggr2ec.pkl', 'rb') as f:
        #     biggr2ec = pickle.load(f)
        # with open('biggec2r.pkl', 'rb') as f:
        #     biggec2r = pickle.load(f)
        ## 如果已经有了负责这个ec的reaction 那 penalty就是 会更高
        # print('orig_rxn_ids:',orig_rxn_ids)
        # for orir in orig_rxn_ids:
        #     ec = biggr2ec[orir]
        #     for curr_ec in ec:
        #         havers = biggec2r[curr_ec]
        #         for nrs in havers:
        #             # remove nrs in weight_dict
        #             if nrs in weight_dict:
        #                 del weight_dict[nrs]
                    # weight_dict[nrs] = 2.0
            # havers = biggec2r[ec]
            # for nrs in havers:
            #     universal_scoredict[nrs] = 0.0
                
        min_old = min(weight_dict.values())
        max_old = max(weight_dict.values())

        pfba_expr = symengine.RealDouble(0)
        for rxn in universal.reactions:
            if rxn.id in orig_rxn_ids:
                coefficientDict[rxn.forward_variable] = 0.0
                coefficientDict[rxn.reverse_variable] = 0.0
            else:
                if rxn.id in weight_dict:
                    # print('rxn.id:',rxn.id,'weight_dict[rxn.id]:',weight_dict[rxn.id],rxn.id in rm_rxn)
                    # coefficientDict[rxn.forward_variable] = weight_dict[rxn.id]
                    # coefficientDict[rxn.reverse_variable] = weight_dict[rxn.id]
                    # coefficientDict[rxn.forward_variable] = 90 * (1 - (weight_dict[rxn.id] - min_old) / (max_old - min_old))
                    # coefficientDict[rxn.reverse_variable] = 90 * (1 - (weight_dict[rxn.id] - min_old) / (max_old - min_old))
                    # coefficientDict[rxn.forward_variable] = 1 * 1/weight_dict[rxn.id]
                    # coefficientDict[rxn.reverse_variable] = 1 * 1/weight_dict[rxn.id]
                    coefficientDict[rxn.forward_variable] = 1 - weight_dict[rxn.id]
                    coefficientDict[rxn.reverse_variable] = 1 - weight_dict[rxn.id]
                else:
                    coefficientDict[rxn.forward_variable] = 1
                    coefficientDict[rxn.reverse_variable] = 1
        # sss =  
        # print('score:',sss)
        # effrs = set()
        # for orir in orig_rxn_ids:
        #     ec = biggr2ec[orir]
        #     for curr_ec in ec:
        #         havers = biggec2r[curr_ec]
        #         for nrs in havers:
        #             # remove nrs in weight_dict
        #             if nrs in weight_dict:
        #                 effrs.add(nrs)
        #                 # try:
        #                 #     coefficientDict[universal.reactions.get_by_id(nrs).forward_variable] += 2
        #                 #     coefficientDict[universal.reactions.get_by_id(nrs).reverse_variable] += 2
        #                 # except:
        #                 #     continue
        # # print('effrs:',len(effrs),effrs)
        # for effr in effrs:
        #     if effr in rm_rxn:
        #         print('effr:',effr)
        #     try:
        #         coefficientDict[universal.reactions.get_by_id(effr).forward_variable] += 10
        #         coefficientDict[universal.reactions.get_by_id(effr).reverse_variable] += 10
        #     except:
        #         continue 
  

        # Create objective, based on pFBA
        print('#create objective, based on pFBA')
        universal.objective = 0
        universal.solver.update()
        universal.objective = universal.problem.Objective(symengine.RealDouble(0), direction='min', sloppy=True)
        universal.objective.set_linear_coefficients(coefficientDict)
        
        # stdout.write('\r[----------------------------------       ]')
        # stdout.flush()

        # Run FBA and identify reactions from universal that are now active
        print('#run FBA and identify reactions from universal that are now active')
        solution = universal.optimize()
        print('# run FBA solution:',solution)
    new_rxn_ids = set([rxn.id for rxn in reaction_bag.reactions if abs(solution.fluxes[rxn.id]) > 1e-6]).difference(orig_rxn_ids)
    # new_rxn_ids = set([rxn.id for rxn in reaction_bag.reactions if abs(solution.fluxes[rxn.id]) > 1e-10]).difference(orig_rxn_ids)
    # stdout.write('\r[-----------------------------------------]\n')
    warnings.filterwarnings('default')

    return(new_rxn_ids) ,coefficientDict


# Add new reactions to model
def _gapfill_model(model, universal, new_rxn_ids, obj, step):
    '''Adds new reactions to model by getting reactions and metabolites to be added to the model, creates gapfilled model, 
    and identifies extracellular metabolites that still need exchanges '''
    # Get reactions and metabolites to be added to the model
    print('#Get reactions and metabolites to be added to the model')
    new_rxns = []
    tmpuni = deepcopy(universal)
    if step == 1: 
        new_rxns.append(tmpuni.reactions.get_by_id(obj).copy())
    for rxn in new_rxn_ids: 
        if rxn != obj:
            new_rxns.append(tmpuni.reactions.get_by_id(rxn).copy())
    
    # Create gapfilled model 
    print('#Create gapfilled model')
    model.add_reactions(new_rxns)
    model.objective = model.problem.Objective(model.reactions.get_by_id(obj).flux_expression, direction='max')

    # Identify extracellular metabolites still need exchanges
    print('#Identify extracellular metabolites still need exchanges')
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
 
def _set_base_inputs(model, universal):
    print('#Set base inputs')
    tasks = ['EX_cpd00035_e','EX_cpd00051_e','EX_cpd00132_e','EX_cpd00041_e','EX_cpd00084_e','EX_cpd00053_e','EX_cpd00023_e',
    'EX_cpd00033_e','EX_cpd00119_e','EX_cpd00322_e','EX_cpd00107_e','EX_cpd00039_e','EX_cpd00060_e','EX_cpd00066_e','EX_cpd00129_e',
    'EX_cpd00054_e','EX_cpd00161_e','EX_cpd00065_e','EX_cpd00069_e','EX_cpd00156_e','EX_cpd00027_e','EX_cpd00149_e','EX_cpd00030_e',
    'EX_cpd00254_e','EX_cpd00971_e','EX_cpd00063_e','EX_cpd10515_e','EX_cpd00205_e','EX_cpd00099_e']

    biggtasks = []
    for task in tasks:
        try:
            test = model.reactions.get_by_id(task)
            biggtasks.append(task)
        except:
            continue
    tasks = biggtasks

    new_rxns = []
    for exch in tasks: 
        try:
            test = model.reactions.get_by_id(exch)
        except:
            new_rxns.append(universal.reactions.get_by_id(exch).copy())
    model.add_reactions(new_rxns)
    for exch in tasks: model.reactions.get_by_id(exch).bounds = (-1000., -0.01)

    return model

def _add_annotation(model, gram, obj='built'):
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