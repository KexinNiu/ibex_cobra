import pickle
import random
import mergem
import mergem
from mergem import translate, load_model, save_model
import pandas as pd
from funcarve_utils import clean2biggr, _find_reactions,_gapfill_model,_set_base_inputs,_add_annotation,_checkModel,read_cleandf_withscore,clean2seedr,weighted_find_reactions
from copy import deepcopy
import argparse
import cobra
import os

# def get_unreviewed_genes(unreviewed_f='/ibex/user/niuk0a/funcarve/cobra/269799_unreviewed.tsv'):
def get_unreviewed_genes(unreviewed_f='../269799_unreviewed.tsv'):
    unreviewed = pd.read_csv(unreviewed_f, sep='\t')
    unreviewed_genes = set(unreviewed['Gene Names'].str.split(' ').sum())
    return unreviewed_genes

def get_rm_pool(model,unreviewed_genes,type='all'):
    can_remove = []
    genes_rxn = {}
    rxn_genes = {}
    for rxn in model.reactions:
        if rxn.gene_reaction_rule != '':
            if 'or' in rxn.gene_reaction_rule or 'and' in rxn.gene_reaction_rule:
                if 'and' in rxn.gene_reaction_rule:
                    genes = rxn.gene_reaction_rule.split(' and ')
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

                        if gene in genes_rxn:
                            genes_rxn[gene].append(rxn.id)
                        else:
                            genes_rxn[gene] = [rxn.id]
                        
                        if gene not in unreviewed_genes:
                            flag = 1
                    if flag == 0:
                        can_remove.append(rxn)
    return can_remove,genes_rxn,rxn_genes


def get_remove_rxn(seednum,prec,can_remove):
    num_rm_rxn = int(len(can_remove) * float(prec))
    random.seed(int(seednum))
    rm_rxn = random.sample(can_remove, num_rm_rxn)
    rm_rxn_id = [rxn.id for rxn in rm_rxn]
    return rm_rxn,rm_rxn_id


def add_rm_rxn_2uni(universal,rm_rxn):
    for r in rm_rxn:
        if r not in universal.reactions:
            print('r:',r.id,'not in universal')
            universal.add_reactions([r.copy()])
    return universal

def gapfill_ori(model,universal_seed,rm_rxn_id,result):
    metabolic_tasks=[]
    gram_type='negative'
    universal_obj = str(model.objective.expression).split()[0].split('*')[-1]
    min_frac=0.01
    max_frac=0.5
    file_type = 1
    exchange_arg=1
    draft_reactions = set([x.id for x in model.reactions])
    draft_metabolites = set([x.id for x in model.metabolites])
    new_reactions = _find_reactions(model, universal_seed, metabolic_tasks, universal_obj, min_frac, max_frac, 1, file_type)
    print('new_reactions:',new_reactions,len(new_reactions))
    inter = set(new_reactions).intersection(set(rm_rxn_id))
    filled_genre = _gapfill_model(model, universal_seed, new_reactions, universal_obj, 1)
    print('Identifying new metabolism (Step 2 of 2)...')
    filled_genre = _set_base_inputs(filled_genre, universal_seed)
    media_reactions = _find_reactions(filled_genre, universal_seed, metabolic_tasks, universal_obj, min_frac, max_frac, 2, file_type)
    
    final_genre = _gapfill_model(filled_genre, universal_seed, media_reactions, universal_obj, 2)
    final_genre = _add_annotation(final_genre, gram_type)
    for exch in final_genre.exchanges: exch.bounds = (-1000., 1000.)
    allnewid = set(new_reactions).union(set(media_reactions))
    inter = allnewid.intersection(set(rm_rxn_id))
    print('allnewid:',allnewid)
    print('rm_rxn_id:',rm_rxn_id)
    print('inter:',inter)
    print('len(inter):',len(inter))
    result['ori_allgfreactions'].append(len(allnewid))
    result['ori_tp'].append(len(inter))
    result['ori_fn'].append(len([rxn for rxn in rm_rxn_id if rxn not in inter]))
    result['ori_fp'].append(len([rxn for rxn in allnewid if rxn not in inter]))
    result['ori_f1'].append(2*len(inter)/(len(rm_rxn_id) + len(allnewid)))
    print('ori:')
    print(f'tp:{len(inter)},fn:{len([rxn for rxn in rm_rxn_id if rxn not in inter])},fp:{len([rxn for rxn in allnewid if rxn not in inter])}')
    return final_genre,result

def gapfill_weight(model,universal,rm_rxn_id,universal_scoredict,result,method):
    metabolic_tasks=[]
    gram_type='negative'
    universal_obj = str(model.objective.expression).split()[0].split('*')[-1]
    min_frac=0.01
    max_frac=0.5
    file_type = 1
    exchange_arg=1
    filled_genre = deepcopy(model)
    draft_reactions = set([x.id for x in filled_genre.reactions])
    draft_metabolites = set([x.id for x in filled_genre.metabolites])
    new_reactions_ids ,coefficientDict,new_rxn_flux = weighted_find_reactions(filled_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 1, file_type,universal_scoredict,method)
    print('Identifying new metabolism (Step 2 of 2)...')
    filled_genre = _set_base_inputs(filled_genre, universal)
    media_reactions ,coefficientDict,new_rxn_flux = weighted_find_reactions(filled_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 2, file_type,universal_scoredict,method)
    final_genre = _gapfill_model(filled_genre, universal, media_reactions, universal_obj, 2)
    final_genre = _add_annotation(final_genre, gram_type)
    allnewid = set(new_reactions_ids).union(set(media_reactions))
    inter = allnewid.intersection(set(rm_rxn_id))
    result['wgf_allgfreactions'].append(len(allnewid))
    result['wgf_tp'].append(len(inter))
    result['wgf_fn'].append(len([rxn for rxn in rm_rxn_id if rxn not in inter]))
    result['wgf_fp'].append(len([rxn for rxn in allnewid if rxn not in inter]))
    result['method'].append(method)
    result['wgf_f1'].append(2*len(inter)/(len(rm_rxn_id) + len(allnewid)))

    print(f'wgf:{method}')
    print(f'tp:{len(inter)},fn:{len([rxn for rxn in rm_rxn_id if rxn not in inter])},fp:{len([rxn for rxn in allnewid if rxn not in inter])}')
    print('all gapfill rxn weight:',len(allnewid))
    # for rxn in allnewid:
    #     print(rxn,universal_scoredict[rxn],'protein:',universal_scoredict[rxn])
    return final_genre,result,allnewid


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Remove reactions based on gene annotation')
    parser.add_argument('--model', type=str, default='../uniprot/iAF987.xml',help='Input model')
    parser.add_argument('--unreviewedf', type=str, default='../269799_unreviewed.tsv',help='Unreviewed genes')
    # parser.add_argument('--cleanf', type=str, default='/ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_t4_maxsep_df.pkl',help='retrain')
    # parser.add_argument('--cleanf', type=str, default='/ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_maxsep_df.pkl',help='ori')
    # /home/kexin/code/bigg/data/CP000148.1_t4_maxsep_df.pkl
    parser.add_argument('--cleanf', type=str, default='/home/kexin/code/bigg/data/CP000148.1_t4_maxsep_df.pkl',help='retrain')
    parser.add_argument('--seednum', type=int, default=7777,help='Seed number for random sampling')
    parser.add_argument('--prec', type=float,default=0.3, help='Percentage of reactions to remove')
    parser.add_argument('--type', type=str,default='all', help='Type of gene annotation')
    parser.add_argument('--threshold', type=int,default=5, help='Type of gene annotation')
    parser.add_argument('--meth', type=int,default=1, help='Type of gene annotation')
    args = parser.parse_args()

    model = cobra.io.read_sbml_model(args.model)  ## bigg model
    media = model.medium ## bigg media
    unreviewed_genes = get_unreviewed_genes(args.unreviewedf)
    pr2ec,ec2pr,predscore,allec2pr = read_cleandf_withscore(args.cleanf,args.threshold)

    method = int(args.meth)
    seednum = args.seednum
    universal_scoredict,rxns,r2maxecp = clean2biggr(predscore,args.threshold)
    
    can_remove,genes_rxn,rxn_genes = get_rm_pool(model,unreviewed_genes,type=args.type)
    universal = cobra.io.load_json_model("../bigg/universal_model_cobrapy.json")
    # universal_seedf = '/ibex/user/niuk0a/funcarve/cobra/universal_reconori.pickle'
    # universal_seed = pickle.load(open(universal_seedf,'rb'))

    if len(media) != 0:
        media_condition = set(['EX_' + cpd for cpd in media])
        universal_reactions = set([x.id for x in universal.reactions])
        for rxn in universal_reactions:
            if rxn.startswith('EX_') == True:
                universal.reactions.get_by_id(rxn).bounds = (0, 1000.0)
            if rxn in media_condition:
                universal.reactions.get_by_id(rxn).bounds = (-1000.0, 10000)

    result = {
        'prec':[],
        'seednum':[],
        'rm_rxn':[],
        'rm_rxn_id':[],
        'method':[],
        'ori_allgfreactions':[],
        'wgf_allgfreactions':[],
        'ori_tp':[],
        'ori_fn':[],
        'ori_fp':[],
        'ori_f1':[],
        'wgf_tp':[],
        'wgf_fn':[],
        'wgf_fp':[],
        'wgf_f1':[]

    }

    prec = args.prec
    tmpmodel = deepcopy(model)
    rm_rxn,rm_rxn_id = get_remove_rxn(args.seednum,prec,can_remove)
    
    for rxnid in rm_rxn_id:
        tmpmodel.reactions.get_by_id(rxnid).remove_from_model()
    rm_model,rm_metabolits = cobra.manipulation.prune_unused_metabolites(tmpmodel)
    ## add rm reaction to universal model
    universal = add_rm_rxn_2uni(universal,rm_rxn)
    cobra.io.write_sbml_model(model, args.model.split('.xml')[0] + '_rm_' + str(prec) + '.sbml')

    result['prec'].append(prec)
    result['seednum'].append(args.seednum)
    result['rm_rxn'].append(len(rm_rxn_id))
    result['rm_rxn_id'].append(rm_rxn_id)
    orirmmodel = deepcopy(rm_model)
    wfgrmmodel = deepcopy(rm_model)
    # ori_gf_model,result = gapfill_ori(orirmmodel,universal,rm_rxn_id,result)
    # print(f'>>Done with ori gapfill {prec}',flush=True)
    w_gf_model,result,allnewrxnid = gapfill_weight(wfgrmmodel,universal,rm_rxn_id,universal_scoredict,result,method)
    print(f'>>Done with weighted gapfill {prec}',flush=True)
    with open('biggr2ec.pkl', 'rb') as f:
        biggr2ec = pickle.load(f)
    print('threshold:',5)
for rxn in allnewrxnid:
    if rxn in biggr2ec:
        ec = biggr2ec[rxn]
        for e in ec:
            e = e.split(':')[1]
            if e in ec2pr.keys():
                pr = ec2pr[e]
                for p in pr.keys():
                    score = predscore[p][e]
                    print(rxn,e,p,score,sep='|')
    else:
        print('nothing:',rxn)
    # cobra.io.write_sbml_model(ori_gf_model, args.model.split('.xml')[0] + '_ori_gf_' + str(prec) + '.sbml')
    # cobra.io.write_sbml_model(w_gf_model, args.model.split('.xml')[0] + '_w_gf_' + str(prec) + '.sbml')
    # break   
    print('end of prec:',prec,flush=True)

# save result
of = open('../data/result/01result_retraindf_t5.csv','a+') ## other all t8
# of = open('/ibex/user/niuk0a/funcarve/cobra/remove_rxn_exp/01result_ori.csv','a+')
print(f'{prec}\t{seednum}\t{len(rm_rxn_id)}\t{rm_rxn_id}\t{result["method"]}\t{result["ori_allgfreactions"]}\t{result["wgf_allgfreactions"]}\t{result["ori_tp"]}\t{result["ori_fn"]}\t{result["ori_fp"]}\t{result["ori_f1"]}\t{result["wgf_tp"]}\t{result["wgf_fn"]}\t{result["wgf_fp"]}\t{result["wgf_f1"]}',file=of)
of.close()


#  the model will minimize overall flux (and thus costly high enzyme turn-over), but still maximize for a given objective function
# 0.005,weight_dict[rxn.id],100
# {'prec': [0.3], 'seednum': [7777], 'rm_rxn': [118], 'wgf_allgfreactions': [1], 'wgf_tp': [0], 'wgf_fn': [118], 'wgf_fp': [1]}
# 0.005,1 - weight_dict[rxn.id],100
#  {'prec': [0.3], 'seednum': [7777], 'rm_rxn': [118],'wgf_allgfreactions': [148], 'wgf_tp': [16], 'wgf_fn': [102], 'wgf_fp': [132]}
# 0.005,1/1+weight_dict[rxn.id],100
#  {'prec': [0.3], 'seednum': [7777], 'rm_rxn': [118], 'wgf_allgfreactions': [101],'wgf_tp': [3], 'wgf_fn': [115], 'wgf_fp': [98]}
# 0.005,weight_dict[rxn.id],1
# {'prec': [0.3], 'seednum': [7777], 'rm_rxn': [118],wgf_allgfreactions': [12],  'wgf_tp': [1], 'wgf_fn': [117], 'wgf_fp': [11]}
# 0.005,1/1+weight_dict[rxn.id],1
# {'prec': [0.3], 'seednum': [7777], 'rm_rxn': [118], 'wgf_allgfreactions': [12], 'wgf_tp': [1], 'wgf_fn': [117], 'wgf_fp': [11]}
# 0.005,1 - weight_dict[rxn.id],1
#  {'prec': [0.3], 'seednum': [7777], 'rm_rxn': [118],'wgf_allgfreactions': [745],'wgf_tp': [17], 'wgf_fn': [101], 'wgf_fp': [728]}



