import os
import pandas as pd
import json
# from read import read_json
# import json
# import pandas as pd


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
# emptytable = pd.DataFrame(columns=['taxo_id','total_metabolites','total_reactions','total_genes','total_compartments','metabolic_coverage','unconserved_metabolites','consistency','annotation_met','annotation_rxn','annotation_gene','annotation_sbo','total_score'])
# def make_table(folder,svtable,name:str,gapfillflag:str):
  
#     taxoid =[]
#     total_metabolites = []
#     total_reactions = []
#     total_genes = []
#     total_compartments = []
#     metabolic_coverage = []
#     unconserved_metabolites = []
#     consistency = []
#     annotation_met = []
#     annotation_rxn = []
#     annotation_gene = []
#     annotation_sbo = []
#     total_score = []
#     newtotal = []
#     con_stoi=[]
#     con_mass=[]
#     con_charge=[]
#     con_met=[]
#     con_unbound=[]
#     con_unbound_n=[]
    
    
#     for root, dirs, files in os.walk(folder):
#         for file in files:
#             # taxo_id = file.split('_')[1]
#             # gapfill = file.split('_')[0]
#             # print('taxo_id:',taxo_id,'gapfill:',gapfill)
#             if file.endswith(".json") and file.split('_')[0] == gapfillflag:
#                 f = os.path.join(root, file)
#                 # KEEG OUT FILE NAME = '267748.5.PATRIC.ori_ngf_267748.sbml'
#                 taxo_id = file.split('_')[1].replace('.json','')
#                 # gapfill = file.split('_')[0]
#                 # print('taxo_id:',taxo_id,'gapfill:',gapfill)
#                 infodict =  read_json(f)

#                 taxoid.append(taxo_id)
#                 total_metabolites.append(infodict['total_metabolites'])
#                 total_reactions.append(infodict['total_reactions'])
#                 total_genes.append(infodict['total_genes'])
#                 total_compartments.append(infodict['total_compartments'])
#                 metabolic_coverage.append(infodict['metabolic_coverage'])
#                 unconserved_metabolites.append(infodict['unconserved_metabolites'])
#                 consistency.append(infodict['consistency'])
#                 annotation_met.append(infodict['annotation_met'])
#                 annotation_rxn.append(infodict['annotation_rxn'])
#                 annotation_gene.append(infodict['annotation_gene'])
#                 annotation_sbo.append(infodict['annotation_sbo'])
#                 total_score.append(infodict['total_score'])
#                 newtotal.append(infodict['newtotal'])
#                 con_stoi.append(infodict['con_stoi'])
#                 con_mass.append(infodict['con_mass'])
#                 con_charge.append(infodict['con_charge'])
#                 con_met.append(infodict['con_met'])
#                 con_unbound.append(infodict['con_unbound'])
#                 con_unbound_n.append(infodict['con_unbound_n'])

    
#     table = pd.DataFrame({'taxo_id':taxoid,
#                           f'{name}_total_metabolites':total_metabolites,
#                           f'{name}_total_reactions':total_reactions,
#                           f'{name}_total_genes':total_genes,
#                           f'{name}_total_compartments':total_compartments,
#                           f'{name}_metabolic_coverage':metabolic_coverage,
#                           f'{name}_unconserved_metabolites':unconserved_metabolites,
#                           f'{name}_consistency':consistency,
#                           f'{name}_annotation_met':annotation_met,
#                           f'{name}_annotation_rxn':annotation_rxn,
#                           f'{name}_annotation_gene':annotation_gene,
#                           f'{name}_annotation_sbo':annotation_sbo,
#                           f'{name}_total_score':total_score,
#                           f'{name}_newtotal':newtotal,
#                             f'{name}_con_stoi':con_stoi,
#                             f'{name}_con_mass':con_mass,
#                             f'{name}_con_charge':con_charge,
#                             f'{name}_con_met':con_met,
#                             f'{name}_con_unbound':con_unbound,
#                             f'{name}_con_unbound_n':con_unbound_n
#                           })   
                                      
#     #sv table
#     table.to_csv(svtable+f'_{gapfillflag}.csv',index=False)    
#     return table


# folder ='/ibex/user/niuk0a/funcarve/bvbrc/memote/kegg'
# svtable = '/ibex/user/niuk0a/funcarve/bvbrc/memote/kegg_n_stat'

# table = make_table(folder,svtable,'keggout','ngf')
# print('table:',table.shape)
# table = make_table(folder,svtable,'keggout','gf')
# print('table:',table.shape)


# folder ='/ibex/user/niuk0a/funcarve/bvbrc/memote/deepec'
# svtable = '/ibex/user/niuk0a/funcarve/bvbrc/memote/deepec_n_stat'

# table = make_table(folder,svtable,'deepec','ngf')
# print('table:',table.shape)
# table = make_table(folder,svtable,'deepec','gf')
# print('table:',table.shape)

# folder = '/ibex/user/niuk0a/funcarve/bvbrc/memote/clean'
# svtable = '/ibex/user/niuk0a/funcarve/bvbrc/memote/clean_n_stat'
# table = make_table(folder,svtable,'clean','ngf')
# print('table:',table.shape)
# table = make_table(folder,svtable,'clean','gf')
# print('table:',table.shape)

#################################################################################################################
def make_table_reclean(folder,svtable,name:str,gapfillflag:str):
  
    # taxoid =[]
    total_metabolites = []
    total_reactions = []
    total_genes = []
    total_compartments = []
    metabolic_coverage = []
    unconserved_metabolites = []
    consistency = []
    annotation_met = []
    annotation_rxn = []
    annotation_gene = []
    annotation_sbo = []
    total_score = []
    newtotal = []
    con_stoi=[]
    con_mass=[]
    con_charge=[]
    con_met=[]
    con_unbound=[]
    con_unbound_n=[]
    modelnames=[]
    evaltypes=[]
    rewards=[]
    thresholds=[]
    iterations=[]
    
    for root, dirs, files in os.walk(folder):
        for file in files:
            # taxo_id = file.split('_')[1]
            # gapfill = file.split('_')[0]
            # print('taxo_id:',taxo_id,'gapfill:',gapfill)
            if file.endswith(".json"):
                f = os.path.join(root, file)
                print('file:',f)
                # /ibex/user/niuk0a/funcarve/cobra/data/result/memotes/AL009126.3_t4iYO844_allec_R1_T6I1.json
                names=file.split('_')[-4:]
                modelname = names[0][2:]
                evaltype = names[1]
                reward = names[2].replace('R','')
                threshold = names[3].split('I')[0].replace('T','')
                iteration= names[3].split('I')[1].replace('.json','')

                modelnames.append(modelname)
                evaltypes.append(evaltype)
                rewards.append(reward)
                thresholds.append(threshold)
                iterations.append(iteration)
                # print(f'modelname:{modelname},eval:{evaltype},reward:{reward},thres:{threshold},inter:{iteration}')


                infodict =  read_json(f)
                # print(infodict)
                # taxoid.append(taxo_id)

                total_metabolites.append(infodict['total_metabolites'])
                total_reactions.append(infodict['total_reactions'])
                total_genes.append(infodict['total_genes'])
                total_compartments.append(infodict['total_compartments'])
                metabolic_coverage.append(infodict['metabolic_coverage'])
                unconserved_metabolites.append(infodict['unconserved_metabolites'])
                consistency.append(infodict['consistency'])
                annotation_met.append(infodict['annotation_met'])
                annotation_rxn.append(infodict['annotation_rxn'])
                annotation_gene.append(infodict['annotation_gene'])
                annotation_sbo.append(infodict['annotation_sbo'])
                total_score.append(infodict['total_score'])
                newtotal.append(infodict['newtotal'])
                con_stoi.append(infodict['con_stoi'])
                con_mass.append(infodict['con_mass'])
                con_charge.append(infodict['con_charge'])
                con_met.append(infodict['con_met'])
                con_unbound.append(infodict['con_unbound'])
                con_unbound_n.append(infodict['con_unbound_n'])

            
    # table = pd.DataFrame({'taxo_id':taxoid,
    table = pd.DataFrame({'modelname':modelnames,
                            'evaltype':evaltypes,
                            'reward':rewards,
                            'threshold':thresholds,
                            'iteration':iterations,
                            'total_metabolites':total_metabolites,
                            'total_reactions':total_reactions,
                            'total_genes':total_genes,
                            'total_compartments':total_compartments,
                            'metabolic_coverage':metabolic_coverage,
                            'unconserved_metabolites':unconserved_metabolites,
                            'consistency':consistency,
                            'annotation_met':annotation_met,
                            'annotation_rxn':annotation_rxn,
                            'annotation_gene':annotation_gene,
                            'annotation_sbo':annotation_sbo,
                            'total_score':total_score,
                            'newtotal':newtotal,
                            'con_stoi':con_stoi,
                            'con_mass':con_mass,
                            'con_charge':con_charge,
                            'con_met':con_met,
                            'con_unbound':con_unbound,
                            'con_unbound_n':con_unbound_n
                        })  
                        #   f'{name}_total_metabolites':total_metabolites,
                        #   f'{name}_total_reactions':total_reactions,
                        #   f'{name}_total_genes':total_genes,
                        #   f'{name}_total_compartments':total_compartments,
                        #   f'{name}_metabolic_coverage':metabolic_coverage,
                        #   f'{name}_unconserved_metabolites':unconserved_metabolites,
                        #   f'{name}_consistency':consistency,
                        #   f'{name}_annotation_met':annotation_met,
                        #   f'{name}_annotation_rxn':annotation_rxn,
                        #   f'{name}_annotation_gene':annotation_gene,
                        #   f'{name}_annotation_sbo':annotation_sbo,
                        #   f'{name}_total_score':total_score,
                        #   f'{name}_newtotal':newtotal,
                        #     f'{name}_con_stoi':con_stoi,
                        #     f'{name}_con_mass':con_mass,
                        #     f'{name}_con_charge':con_charge,
                        #     f'{name}_con_met':con_met,
                        #     f'{name}_con_unbound':con_unbound,
                        #     f'{name}_con_unbound_n':con_unbound_n
                        #   })   
                                      
    #sv table
    table.to_csv(svtable+'.csv',index=False)    
    return table


# folder = '../data/result/memotes/'
# svtable = '../data/result/eval_memote_7v4'
# # table = make_table_reclean(folder,svtable,'cleanr','ngf')
# table = make_table_reclean(folder,svtable,'cleanr','')
# print('table:',table.shape)


def make_table_up(folder,svtable,name:str,gapfillflag:str):
  
    # taxoid =[]
    total_metabolites = []
    total_reactions = []
    total_genes = []
    total_compartments = []
    metabolic_coverage = []
    unconserved_metabolites = []
    consistency = []
    annotation_met = []
    annotation_rxn = []
    annotation_gene = []
    annotation_sbo = []
    total_score = []
    newtotal = []
    con_stoi=[]
    con_mass=[]
    con_charge=[]
    con_met=[]
    con_unbound=[]
    con_unbound_n=[]
    modelnames=[]
    evaltypes=[]
    rewards=[]
    thresholds=[]
    iterations=[]
    
    for root, dirs, files in os.walk(folder):
        for file in files:
            if 'UP' not in file:
                continue
            # taxo_id = file.split('_')[1]
            # gapfill = file.split('_')[0]
            # print('taxo_id:',taxo_id,'gapfill:',gapfill)
            if file.endswith(".json"):
                f = os.path.join(root, file)
                print('file:',f)
                # /ibex/user/niuk0a/funcarve/cobra/data/result/memotes/AL009126.3_t4iYO844_allec_R1_T6I1.json
                #/ibex/user/niuk0a/funcarve/cobra/data/result/memotes/UP000000640_196162I4.json
                names=file.split('_')[-4:]
                modelname = names[0]
                taxoid = names[1][:-2]
                
                iteration= names[1].split('I')[1].replace('.json','')
                modelnames.append(modelname)
                evaltypes.append('differ')
                
                iterations.append(iteration)
                # print(f'modelname:{modelname},eval:{evaltype},reward:{reward},thres:{threshold},inter:{iteration}')


                infodict =  read_json(f)
                # print(infodict)
                # taxoid.append(taxo_id)

                total_metabolites.append(infodict['total_metabolites'])
                total_reactions.append(infodict['total_reactions'])
                total_genes.append(infodict['total_genes'])
                total_compartments.append(infodict['total_compartments'])
                metabolic_coverage.append(infodict['metabolic_coverage'])
                unconserved_metabolites.append(infodict['unconserved_metabolites'])
                consistency.append(infodict['consistency'])
                annotation_met.append(infodict['annotation_met'])
                annotation_rxn.append(infodict['annotation_rxn'])
                annotation_gene.append(infodict['annotation_gene'])
                annotation_sbo.append(infodict['annotation_sbo'])
                total_score.append(infodict['total_score'])
                newtotal.append(infodict['newtotal'])
                con_stoi.append(infodict['con_stoi'])
                con_mass.append(infodict['con_mass'])
                con_charge.append(infodict['con_charge'])
                con_met.append(infodict['con_met'])
                con_unbound.append(infodict['con_unbound'])
                con_unbound_n.append(infodict['con_unbound_n'])

            
    # table = pd.DataFrame({'taxo_id':taxoid,
    table = pd.DataFrame({'modelname':modelnames,
                            'evaltype':evaltypes,
                            'iteration':iterations,
                            'total_metabolites':total_metabolites,
                            'total_reactions':total_reactions,
                            'total_genes':total_genes,
                            'total_compartments':total_compartments,
                            'metabolic_coverage':metabolic_coverage,
                            'unconserved_metabolites':unconserved_metabolites,
                            'consistency':consistency,
                            'annotation_met':annotation_met,
                            'annotation_rxn':annotation_rxn,
                            'annotation_gene':annotation_gene,
                            'annotation_sbo':annotation_sbo,
                            'total_score':total_score,
                            'newtotal':newtotal,
                            'con_stoi':con_stoi,
                            'con_mass':con_mass,
                            'con_charge':con_charge,
                            'con_met':con_met,
                            'con_unbound':con_unbound,
                            'con_unbound_n':con_unbound_n
                        })  
                        #   f'{name}_total_metabolites':total_metabolites,
                        #   f'{name}_total_reactions':total_reactions,
                        #   f'{name}_total_genes':total_genes,
                        #   f'{name}_total_compartments':total_compartments,
                        #   f'{name}_metabolic_coverage':metabolic_coverage,
                        #   f'{name}_unconserved_metabolites':unconserved_metabolites,
                        #   f'{name}_consistency':consistency,
                        #   f'{name}_annotation_met':annotation_met,
                        #   f'{name}_annotation_rxn':annotation_rxn,
                        #   f'{name}_annotation_gene':annotation_gene,
                        #   f'{name}_annotation_sbo':annotation_sbo,
                        #   f'{name}_total_score':total_score,
                        #   f'{name}_newtotal':newtotal,
                        #     f'{name}_con_stoi':con_stoi,
                        #     f'{name}_con_mass':con_mass,
                        #     f'{name}_con_charge':con_charge,
                        #     f'{name}_con_met':con_met,
                        #     f'{name}_con_unbound':con_unbound,
                        #     f'{name}_con_unbound_n':con_unbound_n
                        #   })   
                                      
    #sv table
    table.to_csv(svtable+'.csv',index=False)    
    return table


folder = '../data/result/memotes/'
svtable = '../data/result/eval_memote_t2up'
# table = make_table_reclean(folder,svtable,'cleanr','ngf')
table = make_table_up(folder,svtable,'cleanr','')
print('table:',table.shape)



