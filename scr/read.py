import json
import pandas as pd


def read_json(f):
    with open(f) as json_file:
        data = json.load(json_file)
    infodict = {}
    infodict['total_score'] = data['score']['total_score']
    infodict['total_metabolites'] = data['tests']['test_metabolites_presence']['message'].split(' ')[0]
    infodict['total_reactions'] = data['tests']['test_reactions_presence']['message'].split(' ')[0]
    genes_messages =data['tests']['test_genes_presence']['message'].split(' ')
    infodict['total_genes'] = genes_messages[0]
    infodict['total_compartments'] = data['tests']['test_compartments_presence']['message'].split(' ')[3]
    try:
        infodict['metabolic_coverage'] = data['tests']['test_metabolic_coverage']['message'].split(' ')[6][:-1]
    except:
        infodict['metabolic_coverage'] = data['tests']['test_metabolic_coverage']['metric']
    infodict['unconserved_metabolites'] = data['tests']['test_unconserved_metabolites']['message'].split(' ')[3]
    infodict['consistency'] = data['score']['sections'][0]['score']
    
   
    stoi = data['tests']["test_stoichiometric_consistency"]['data']
    if stoi:
        infodict['con_stoi'] =1
    else:
        infodict['con_stoi'] =0
    try:    
        infodict['con_mass'] = 100-float(data['tests']["test_reaction_mass_balance"]['message'].split(' ')[4][1:-2])
    except:
        infodict['con_mass'] = 0
    try:    
        infodict['con_charge'] = 100-float(data['tests']["test_reaction_charge_balance"]['message'].split(' ')[4][1:-2])
    except:
        infodict['con_charge'] = 0
    try:
        infodict['con_met'] = 100-float(data['tests']["test_find_disconnected"]['message'].split(' ')[4][1:-2])
    except:
        infodict['con_met'] = 0
    try:
        infodict['con_unbound'] = 100-float(data['tests']["test_find_reactions_unbounded_flux_default_condition"]['message'].split(' of ')[1][0:-1])
    except:
        infodict['con_unbound'] = 0
    try:
        infodict['con_unbound_n'] = data['tests']["test_find_reactions_unbounded_flux_default_condition"]['message'].split('in total ')[1].split('\nreactions')[0]
    except:
        infodict['con_unbound_n'] = 0
    infodict['annotation_met'] = data['score']['sections'][1]['score']
    infodict['annotation_rxn'] = data['score']['sections'][2]['score']
    infodict['annotation_gene'] = data['score']['sections'][3]['score']
    infodict['annotation_sbo'] = data['score']['sections'][4]['score']
    infodict['newtotal'] = (float(infodict['consistency'])*3 + float(infodict['annotation_met']) +  float(infodict['annotation_rxn']) +2* float(infodict['annotation_sbo'])) /7
    return infodict

# f = '/ibex/user/niuk0a/funcarve/bvbrc/memote/mgen_dgm_euc7.json'
# infodict = read_json(f)
# print(infodict)
# print(infodict.keys())  
