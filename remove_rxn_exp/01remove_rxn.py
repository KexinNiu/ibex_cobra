import pickle
import random
# import mergem
# from mergem import translate, load_model, save_model
import pandas as pd
from rm_funcarve_utils import clean2biggr, find_reactions,gapfill_model,set_base_inputs,add_annotation,checkModel,read_cleandf_withscore,clean2seedr,weighted_find_reactions
from Bio import SeqIO
from copy import deepcopy
import argparse
import cobra
import os
import pandas as pd



def parse_genebank(f):
    recs = [rec for rec in SeqIO.parse(f, "genbank")]
    rec = recs[0]
    feats = [feat for feat in rec.features if feat.type == "CDS"]
    lt2ec={}
    lt2oldlt={}
    oldlt2lt={}
    for feat in feats:
        dd = feat.qualifiers
        '''
        dd = {'locus_tag': ['TM_RS00005'], 'old_locus_tag': ['TM0005', 'TM_0005'], 'EC_number': ['3.6.4.12'], 'inference': ['COORDINATES: similar to AA sequence:RefSeq:WP_012310830.1'], 'GO_function': ['GO:0003678 - DNA helicase activity [Evidence IEA]'], 'GO_process': ['GO:0006281 - DNA repair [Evidence IEA]'], 'note': ['Derived by automated computational analysis using gene prediction method: Protein Homology.'], 'codon_start': ['1'], 'transl_table': ['11'], 'product': ['IGHMBP2 family helicase'], 'protein_id': ['WP_010865024.1'], 'db_xref': ['GI:499163180'], 'translation': ['MTVQQFIKKLVRLVELERNAEINAMLDEMKRLSGEEREKKGRAVLGLTGKFIGEELGYFLVRFGRRKKIDTEIGVGDLVLISKGNPLKSDYTGTVVEKGERFITVAVDRLPSWKLKNVRIDLFASDITFRRQIENLMTLSSEGKKALEFLLGKRKPEESFEEEFTPFDEGLNESQREAVSLALGSSDFFLIHGPFGTGKTRTLVEYIRQEVARGKKILVTAESNLAVDNLVERLWGKVSLVRIGHPSRVSSHLKESTLAHQIETSSEYEKVKKMKEELAKLIKKRDSFTKPSPQWRRGLSDKKILEYAEKNWSARGVSKEKIKEMAEWIKLNSQIQDIRDLIERKEEIIASRIVREAQVVLSTNSSAALEILSGIVFDVVVVDEASQATIPSILIPISKGKKFVLAGDHKQLPPTILSEDAKDLSRTLFEELITRYPEKSSLLDTQYRMNELLMEFPSEEFYDGKLKAAEKVRNITLFDLGVEIPNFGKFWDVVLSPKNVLVFIDTKNRSDRFERQRKDSPSRENPLEAQIVKEVVEKLLSMGVKEDWIGIITPYDDQVNLIRELIEAKVEVHSVDGFQGREKEVIIISFVRSNKNGEIGFLEDLRRLNVSLTRAKRKLIATGDSSTLSVHPTYRRFVEFVKKKGTYVIF']}
        '''
        locus_tag = dd['locus_tag']
        try:
            ec_number = dd['EC_number']
            lt2ec[locus_tag[0]] = ec_number
        except:
            pass
        try:
            old_locus_tag = dd['old_locus_tag']
            lt2oldlt[locus_tag[0]] = old_locus_tag
            if len(old_locus_tag) >= 1:
                for old in old_locus_tag:
                    oldlt2lt[old] = locus_tag
        except:
            pass
        
    return lt2ec, lt2oldlt, oldlt2lt   
 
# def get_unreviewed_genes(unreviewed_f='/ibex/user/niuk0a/funcarve/cobra/269799_unreviewed.tsv'):
def get_reviewed_genes(reviewed_f):
    reviewed = pd.read_csv(reviewed_f, sep='\t')
    reviewed = reviewed[reviewed['Reviewed'] == 'reviewed']
    ## drop nan
    reviewed = reviewed.dropna(subset=['Gene Names'])
    reviewed_genes = set(reviewed['Gene Names'].str.split(' ').sum())
    return reviewed_genes

def get_rm_pool(model,reviewed_genes,type='all'):
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
                    
                        if  gene in reviewed_genes:
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
                        if gene in reviewed_genes:
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
                        
                        if gene not in reviewed_genes:
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

    new_reactions = find_reactions(model, universal_seed, metabolic_tasks, universal_obj, min_frac, max_frac, 1, file_type)
    print('new_reactions:',new_reactions,len(new_reactions))
    inter = set(new_reactions).intersection(set(rm_rxn_id))
    filled_genre = gapfill_model(model, universal_seed, new_reactions, universal_obj, 1)
    print('Identifying new metabolism (Step 2 of 2)...')
    filled_genre = set_base_inputs(filled_genre, universal_seed)
    media_reactions = find_reactions(filled_genre, universal_seed, metabolic_tasks, universal_obj, min_frac, max_frac, 2, file_type)
    
    final_genre = gapfill_model(filled_genre, universal_seed, media_reactions, universal_obj, 2)
    final_genre = add_annotation(final_genre, gram_type)
    for exch in final_genre.exchanges: exch.bounds = (-1000., 1000.)
    allnewid = set(new_reactions).union(set(media_reactions))
    inter = allnewid.intersection(set(rm_rxn_id))
    print('ORI:')
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
    if method ==1:
        upper = 15
        lower = 5
        maxw =100.0
        minw =0.0
    elif method ==2:
        upper = 15
        lower = 6
        maxw = 100.0
        minw = 0.0
    elif method ==3:
        upper = 20
        lower = 5
        maxw = 100.0
        minw = 0.0
    elif method ==4:
        upper = 20
        lower = 4
        maxw = 100.0
        minw = 0.0
    elif method ==5:
        upper = 20
        lower = 4
        maxw = 100.0
        minw = 0.05
    elif method ==6:
        upper = 15
        lower = 6
        maxw = 100.0
        minw = 0.05

    file_type = 1
    exchange_arg=1
    filled_genre = deepcopy(model)
    draft_reactions = set([x.id for x in filled_genre.reactions])
    draft_metabolites = set([x.id for x in filled_genre.metabolites])
    new_reactions_ids ,coefficientDict,new_rxn_flux = weighted_find_reactions(filled_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 1, file_type,universal_scoredict,upper,lower,maxw,minw)
    print('Identifying new metabolism (Step 2 of 2)...')
    filled_genre = set_base_inputs(filled_genre, universal)
    media_reactions ,coefficientDict,new_rxn_flux = weighted_find_reactions(filled_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 2, file_type,universal_scoredict,upper,lower,maxw,minw)
    final_genre = gapfill_model(filled_genre, universal, media_reactions, universal_obj, 2)
    final_genre = add_annotation(final_genre, gram_type)
    allnewid = set(new_reactions_ids).union(set(media_reactions))
    inter = allnewid.intersection(set(rm_rxn_id))
    print('WGF:')
    print('rm_rxn_id:',len(rm_rxn_id))
    result['wgf_allgfreactions'].append(len(allnewid))
    result['wgf_tp'].append(len(inter))
    result['wgf_fn'].append(len([rxn for rxn in rm_rxn_id if rxn not in inter]))
    result['wgf_fp'].append(len([rxn for rxn in allnewid if rxn not in inter]))
    result['method'].append(method)
    result['wgf_f1'].append(2*len(inter)/(len(rm_rxn_id) + len(allnewid)))

    print(f'wgf:{method}')
    print(f'tp:{len(inter)},fn:{len([rxn for rxn in rm_rxn_id if rxn not in inter])},fp:{len([rxn for rxn in allnewid if rxn not in inter])}')
    print('all gapfill rxn weight:',len(allnewid))
    return final_genre,result,allnewid


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Remove reactions based on gene annotation')
    parser.add_argument('--model', type=str, default='../uniprot/iAF987.xml',help='Input model')
    parser.add_argument('--unreviewedf', type=str, default='../269799_unreviewed.tsv',help='Unreviewed genes')
    # parser.add_argument('--cleanf', type=str, default='/ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_t4_maxsep_df.pkl',help='retrain')
    parser.add_argument('--cleanf', type=str, default='/home/kexin/code/bigg/data/CP000148.1_t4_maxsep_df.pkl',help='retrain')
    parser.add_argument('--seednum', type=int, default=7777,help='Seed number for random sampling')
    parser.add_argument('--prec', type=float,default=0.3, help='Percentage of reactions to remove')
    parser.add_argument('--type', type=str,default='all', help='Type of gene annotation')
    parser.add_argument('--meth', type=int,default=1, help='Type of gene annotation')
    parser.add_argument('--ori',type=int,default=1, help='Type of gene annotation')
    args = parser.parse_args()
######
    # python 01remove_rxn.py --prec 0.2  --seednum 7777 --meth 1 --cleanf  /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl --unreviewedf /ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_243274.tsv --model /ibex/user/niuk0a/funcarve/cobra/uniprot/iLJ478.xml
    # method = 1
    # seednum = 7777
    # args.cleanf = '/ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl'
    # args.unreviewedf = '/ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_243274.tsv'
    # args.model = '/ibex/user/niuk0a/funcarve/cobra/uniprot/iLJ478.xml'
    # args.prec = 0.2
    # args.meth = 1
    modelname = args.model.split('/')[-1].split('.')[0]
    print('modelname:',modelname)
    ori = args.ori
    print('ori:',ori)

######

    model = cobra.io.read_sbml_model(args.model)  ## bigg model
    media = model.medium ## bigg media
    reviewed_genes = get_reviewed_genes(args.unreviewedf)
    # reviewed_genes need to translate by .gb file
    gbf = args.cleanf.replace('_t1_maxsep_df.pkl','.gb').replace('/ibex/user/niuk0a/CLEAN/app/results/inputs/','/ibex/user/niuk0a/funcarve/cobra/uniprot/')
    print('gbf:',gbf)
    print('ori reviewed_genes:',len(reviewed_genes),list(reviewed_genes)[:20])
    new_reviewed_genes = set()
    if os.path.exists(gbf):
        lt2ec, lt2oldlt, oldlt2lt  = parse_genebank(gbf)
        for gene in reviewed_genes:
            if gene in oldlt2lt:
                # reviewed_genes.remove(gene)
                for item in oldlt2lt[gene]:
                    new_reviewed_genes.add(item)
                    otherold = lt2oldlt[item]
                    for old in otherold:
                        new_reviewed_genes.add(old)

    reviewed_genes = new_reviewed_genes.union(reviewed_genes)
                # reviewed_genes.add(oldlt2lt[gene])

    print('reviewed_genes:',len(reviewed_genes),list(reviewed_genes)[:10])
    pr2ec,ec2pr,predscore,allec2pr = read_cleandf_withscore(args.cleanf,threshold=5)
    # universal_scoredict,rxns = clean2seedr(predscore,threshold=8)
    method = int(args.meth)
    seednum = args.seednum





    universal_scoredict,rxns,r2maxecp = clean2biggr(predscore,threshold=5)
    
    can_remove,genes_rxn,rxn_genes = get_rm_pool(model,reviewed_genes,type=args.type)
    print('can_remove:',len(can_remove))
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
         #,
        # 'ori_blockr':[],
        # 'wgf_blockr':[]
    }
    # precrange = [0.1,0.2,0.3,0.4,0.5]
    prec = args.prec
    # for prec in precrange:
    tmpmodel = deepcopy(model)
    rm_rxn,rm_rxn_id = get_remove_rxn(args.seednum,prec,can_remove)
    print('rm_rxn:',rm_rxn)
    print('rm_rxn_id:',rm_rxn_id)
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
    # ori_gf_model,result = gapfill_ori(model,universal_seed,rm_rxn_id,result)
    print('rm_rxn_id:',len(rm_rxn_id))
    if ori == 1:
        ori_gf_model,result = gapfill_ori(orirmmodel,universal,rm_rxn_id,result)
        print(f'>>Done with ori gapfill {prec}',flush=True)
        print('rm_rxn_id:',len(rm_rxn_id))
    else:
        w_gf_model,result,allnewrxnid = gapfill_weight(wfgrmmodel,universal,rm_rxn_id,universal_scoredict,result,method)
        print(f'>>Done with weighted gapfill {prec}',flush=True)
 
    print('end of prec:',prec,flush=True)
    print('result:',[(key,len(item)) for key,item in result.items() ])

# save result
# of = open('01result_retraindf_t5.csv','a+') ## other all t8
# print(f'{prec}\t{seednum}\t{len(rm_rxn_id)}\t{rm_rxn_id}\t{result["method"]}\t{result["ori_allgfreactions"]}\t{result["wgf_allgfreactions"]}\t{result["ori_tp"]}\t{result["ori_fn"]}\t{result["ori_fp"]}\t{result["ori_f1"]}\t{result["wgf_tp"]}\t{result["wgf_fn"]}\t{result["wgf_fp"]}\t{result["wgf_f1"]}',file=of)
# # print(f'{prec}\t{seednum}\t{len(rm_rxn_id)}\t{rm_rxn_id}\t{result["method"]}\t{result["ori_allgfreactions"]}\t{result["wgf_allgfreactions"]}\t{result["ori_tp"]}\t{result["ori_fn"]}\t{result["ori_fp"]}\t{result["ori_f1"]}\t{result["wgf_tp"]}\t{result["wgf_fn"]}\t{result["wgf_fp"]}\t{result["wgf_f1"]}',file=of)
# of.close()

if ori == 1:
    oriof = open(f'mmodel7_{modelname}_ori_r.csv','a+') ## other all t8
    print(f'{modelname}\t{prec}\t{seednum}\t{len(rm_rxn_id)}\t{rm_rxn_id}\t{result["method"]}\t{result["ori_allgfreactions"]}\t{result["ori_tp"]}\t{result["ori_fn"]}\t{result["ori_fp"]}\t{result["ori_f1"]}',file=oriof)
    oriof.close()
else:
    wgf_of = open(f'mmodel7_{modelname}_wgf_r.csv','a+') ## other all t8
    print(f'{modelname}\t{prec}\t{seednum}\t{len(rm_rxn_id)}\t{rm_rxn_id}\t{result["method"]}\t{result["wgf_allgfreactions"]}\t{result["wgf_tp"]}\t{result["wgf_fn"]}\t{result["wgf_fp"]}\t{result["wgf_f1"]}',file=wgf_of)
    wgf_of.close()
# print('result:',[(key,item) for key,item in result.items() if key != 'rm_rxn_id'])
# # exit()

# if os.path.exists(f'remove_rxn_result{seednum}.pickle'):
#     with open(f'remove_rxn_result{seednum}.pickle','rb') as f:
#         resultdd = pickle.load(f)
#     for i in range(0,len(result['prec'])):
#         resultdd['prec'].append(result['prec'][i])
#         resultdd['seednum'].append(result['seednum'][i])
#         resultdd['rm_rxn'].append(result['rm_rxn'][i])
#         resultdd['rm_rxn_id'].append(result['rm_rxn_id'][i])
#         # resultdd['ori_allgfreactions'].append(result['ori_allgfreactions'][i])
#         resultdd['wgf_allgfreactions'].append(result['wgf_allgfreactions'][i])
#         # resultdd['ori_tp'].append(result['ori_tp'][i])
#         # resultdd['ori_fn'].append(result['ori_fn'][i])
#         # resultdd['ori_fp'].append(result['ori_fp'][i])
#         # resultdd['ori_f1'].append(result['ori_f1'][i])
#         resultdd['wgf_tp'].append(result['wgf_tp'][i])
#         resultdd['wgf_fn'].append(result['wgf_fn'][i])
#         resultdd['wgf_fp'].append(result['wgf_fp'][i])
#         resultdd['wgf_f1'].append(result['wgf_f1'][i])
#         resultdd['method'].append(result['method'][i])
#     resultdd = pd.DataFrame(resultdd)
#     with open(f'01result.csv','a+') as of:
#         print(f'{prec}\t{seednum}\t{rm_rxn_id}\t{result["ori_allgfreactions"]}\t{result["wgf_allgfreactions"]}\t{result["ori_tp"]}\t{result["ori_fn"]}\t{result["ori_fp"]}\t{result["ori_f1"]}\t{result["wgf_tp"]}\t{result["wgf_fn"]}\t{result["wgf_fp"]}\t{result["wgf_f1"]}',file=of)
# else:
#     result.pop('ori_allgfreactions')
#     result.pop('ori_tp')
#     result.pop('ori_fn')
#     result.pop('ori_fp')
#     result.pop('ori_f1')

#     resultdd = pd.DataFrame(result)
# with open(f'remove_rxn_result{seednum}.pickle','wb') as f:
#     pickle.dump(resultdd,f)
# resultdd.to_csv(f'remove_rxn_result{seednum}.csv',index=False)


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



