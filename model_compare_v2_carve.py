from copy import deepcopy
import warnings
import cobra
from recon import _create_model,_find_reactions, _gapfill_model, _add_annotation, _checkModel, _set_base_inputs, w_find_reactions,print_model_info,read_clean_withscore,clean2biggr
import re
import argparse
import random
parser = argparse.ArgumentParser(description='Precentage to remove reactions')
parser.add_argument('--model','-m', help='model path')
parser.add_argument('--type', '-t',default = 'gf', help='gf or gf weighted')
parser.add_argument('--cleanfile','-cf', help='gf or gf weighted')
parser.add_argument('--prec', default = 0.05, help='percentage of reactions to remove')
parser.add_argument('--seednum', default = 1234, help='seed for random')

args = parser.parse_args()

def detect_obj(m):
    biomass_re = re.compile("biomass", re.IGNORECASE)
    possible_objectives = m.reactions.query(biomass_re)
    if len(m.reactions.query(lambda x: x > 0, "objective_coefficient")):
        return m.objective_coefficient
    # look for reactions with "biomass" in the id or name
    possible_objectives = m.reactions.query(biomass_re)
    if len(possible_objectives) == 0:
        possible_objectives = m.reactions.query(biomass_re, "name")
    
    # # In some cases, a biomass "metabolite" is produced, whose production
    # # should be the objective function.
    # possible_biomass_metabolites = m.metabolites.query(biomass_re)
    # if len(possible_biomass_metabolites) == 0:
    #     possible_biomass_metabolites = m.metabolites.query(biomass_re, "name")

    # if len(possible_biomass_metabolites) > 0:
    #     biomass_met = possible_biomass_metabolites[0]
    #     r = cobra.Reaction("added_biomass_sink")
    #     r.objective_coefficient = 1
    #     r.add_metabolites({biomass_met: -1})
    #     m.add_reaction(r)
    #     print ("autodetected biomass metabolite '%s' for model '%s'" %
    #           (biomass_met.id, m.id))

    elif len(possible_objectives) > 0:
        print("autodetected objective reaction '%s' for model '%s'" %
              (possible_objectives[0].id, m.id))
        m.change_objective(possible_objectives[0])

    else:
        print("no objective found for " + m.id)
    return m.objective.expression

def bigg_obj(m):
    maxval = 0 
    obj = ''
    items = str(m.objective.expression).split(' ')
    # print('items:',items)
    for item in items:
        if '*' in item:
            if abs(float(item.split('*')[0])) > maxval:
                maxval = abs(float(item.split('*')[0]))
                obj = item.split('*')[1]
            if 'BIOMASS' in item or 'biomass' in item:
                print('biomass obj found:',item)
                obj = item.split('*')[1]
                print('obj:',obj)
                return obj
    return obj

def create_model_randomrm(model,universal,prec,seednum,input_id):
    obj = bigg_obj(model)
    print('bigg_obj:',obj)  
    model = deepcopy(model)

    tmpm = deepcopy(model)
    for r in tmpm.reactions:
        try:
            universal.add_reaction(r.copy())
        except:
            pass
    
    random.seed(seednum)
    rxns = list(model.reactions)
    num_rxns = len(rxns)
    num_rm = int(num_rxns * (1-float(prec)))
    sample_rxn = random.sample(rxns,num_rm)
    sample_rxn_id = [r.id for r in sample_rxn]
    print('r1 =',sample_rxn[0].id)
    print('Number of reactions to remove:',int(num_rxns * float(prec)))
    print('removed rxns:',[r.id for r in rxns if r not in sample_rxn])
    rm_model = _create_model(sample_rxn_id,universal,input_id)
    rm_model.objective = obj
    # print('rm_model:',type(rm_model))
    # print('now obj:',rm_model.objective)
    # exp = rm_model.reactions.get_by_id(obj).flux_expression
    # print('exp:',exp)
    # rm_model.objective = rm_model.problem.Objective(rm_model.reactions.get_by_id(obj).flux_expression, direction='max')
    
    # print('rm_model.objnow:',rm_model.objective)
    # # print('model.objnow:',model.objective)
    # # print('split:',str(model.objective.expression).split()[0].split('*')[-1])
    # # try:
    # #     rm_model.objective = obj
    # # except:
    # #     rm_model.objective = rm_r
    # print('rm_model.objnow:',rm_model.objective)
    # rm_model.objective = str(model.objective.expression).split()[0].split('*')[-1]
    return rm_model
    
    
    
def gapfill_draftmodel(draftmodel,universal,gftype,cleanfile=None):
    exchange_arg = 1 # 0 for no exchange, 1 for exchange
    metabolic_tasks=[] # no metabolic tasks
    min_frac = 0.01 # Minimum objective fraction required during gapfilling
    max_frac = 0.5  # Maximum objective fraction required during gapfilling
    file_type = 3 # 3 for sbml draft model
    
    model = deepcopy(draftmodel)    
    print_model_info(model)
    # universal = cobra.io.load_json_model("/ibex/user/niuk0a/funcarve/cobra/bigg/universal_model_cobrapy.json")
    media = model.medium
    
    if len(media) != 0:
        media_condition = set(['EX_' + cpd for cpd in media])
        universal_reactions = set([x.id for x in universal.reactions])
        for rxn in universal_reactions:
            if rxn.startswith('EX_') == True:
                universal.reactions.get_by_id(rxn).bounds = (0, 1000.0)
            if rxn in media_condition:
                universal.reactions.get_by_id(rxn).bounds = (-1000.0, 10000)
    universal_obj = str(model.objective.expression).split()[0].split('*')[-1]           
    
    if gftype == 'gf':
        print('Gapfilling type: gf')
        # universal_obj = str(model.objective.expression).split()[0].split('*')[-1]
        print('Identifying new metabolism...')
        draft_reactions = set([x.id for x in model.reactions])
        draft_metabolites = set([x.id for x in model.metabolites])
        warnings.filterwarnings('ignore')
        new_reactions = _find_reactions(model, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 1, file_type)
        filled_genre = _gapfill_model(model, universal, new_reactions, universal_obj, 1)
        print('Identifying new metabolism (Step 2 of 2)...')
        filled_genre = _set_base_inputs(filled_genre, universal)
        media_reactions = _find_reactions(filled_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 2, file_type)
        final_genre = _gapfill_model(filled_genre, universal, media_reactions, universal_obj, 2)
        final_genre = _add_annotation(filled_genre, universal_obj)
        
        print('GF new reactions:',len(new_reactions),'|', 'GF media reactions:',len(media_reactions))
        print('GF new reactions:',new_reactions,'\n', 'GF media reactions:',media_reactions)
        
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
        ## save the final model
        cobra.io.write_sbml_model(final_genre, modelf.split('.xml')[0] + '_gf_v1.xml')
        
    elif gftype == 'gfw':
        print('Gapfilling type: gfw')
        # universal_obj = str(model.objective.expression).split()[0].split('*')[-1]
        pr2ec , predscore = read_clean_withscore(cleanfile)
        # pr2ec , predscore = read_clean_withscore('/ibex/user/niuk0a/funcarve/cobra/aaseq_CP000148.1_maxsep.csv')
        universal_scoredict,r2maxecp = clean2biggr(predscore) # r2maxecp[rid]=[pr,ec,score]
        print('Identifying new metabolism...')
        draft_reactions = set([x.id for x in model.reactions])
        draft_metabolites = set([x.id for x in model.metabolites])
        warnings.filterwarnings('ignore')
        new_reactions_ids, coefficientDict = w_find_reactions(model, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 1, file_type,universal_scoredict)
        print('new reactions:',len(new_reactions_ids))
        filled_genre = _gapfill_model(model, universal, new_reactions_ids, universal_obj, 1)
        for rid in new_reactions_ids:
            print(f"reaction:{rid} add with {r2maxecp[rid]}")
        print('Identifying new metabolism (Step 2 of 2)...')    
        filled_genre = _set_base_inputs(filled_genre, universal)
        media_reactions,coefficientDict = w_find_reactions(filled_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 2, file_type,universal_scoredict)
        final_genre = _gapfill_model(filled_genre, universal, media_reactions, universal_obj, 2)
        final_genre = _add_annotation(filled_genre, universal_obj)
        print('GFW new reactions:',len(new_reactions_ids),'|', 'GFW media reactions:',len(media_reactions))
        print('GFW new reactions:',new_reactions_ids,'\n', 'GFW media reactions:',media_reactions)
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
        ## save the final model
        cobra.io.write_sbml_model(final_genre, modelf.split('.xml')[0] + '_gfw_v1.xml')
    
    
if __name__ == '__main__':
    universal = cobra.io.load_json_model("/ibex/user/niuk0a/funcarve/cobra/bigg/universal_model_cobrapy.json")
    modelf = args.model
    gftype = args.type
    prec = args.prec
    seednum = args.seednum
    input_id = str(args.model).split('/')[-1].split('.xml')[0] + '_rm' + str(gftype)
    
    
    if gftype == 'gf':
        cleanfile = None
    elif gftype == 'gfw':
        cleanfile = args.cleanfile
        if cleanfile == None:
            print('Please provide the cleanfile')
            exit()
    m = cobra.io.read_sbml_model(modelf)
    model = create_model_randomrm(m,universal,prec,seednum,input_id)
    
    gapfill_draftmodel(model,universal,gftype,cleanfile)
    print('Done!')