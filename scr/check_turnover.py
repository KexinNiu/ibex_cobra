import pickle
import pandas as pd
import csv
import numpy as np
import os
from Bio import SeqIO
import argparse
# from eval_utils import get_metric, get_ancester
from funcarve_utils import differ

def get_ancester(eclist):
    allecs = set()
    for ec in eclist:
        digs =  ec.split('.')
        for i in range(4):
            allecs.add('.'.join(digs[:i+1]))
    return allecs

def get_metric(ec_pred,ec_true):
    ec_pred = get_ancester(ec_pred)
    ec_true = get_ancester(ec_true)
    # print('type ec_pred',type(ec_pred))
    # print('type ec_true',type(ec_true))
    # print('all len',len(ec_pred),len(ec_true))
    union_len = len(ec_pred.union(ec_true))
    true_len = len(ec_true)
    pred_len = len(ec_pred)
    
    acc = len(ec_pred.intersection(ec_true)) / union_len if union_len != 0 else 0
    recall = len(ec_pred.intersection(ec_true)) / true_len if true_len != 0 else 0
    precison = len(ec_pred.intersection(ec_true)) / pred_len if pred_len != 0 else 0
    f1_score = 2 * recall * precison / (recall + precison) if (recall + precison) != 0 else 0
    # print('acc',acc,'recall',recall,'precison',precison,'f1_score',f1_score)
    return acc, recall, precison,f1_score


# orif='test_cp148_t5.0_ORIpredscore_0.pkl'
# with open(orif, 'rb') as file:
#     oridata = pickle.load(file)

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

def read_unif_genebk(f,gbkf):
    print('gbkf func',gbkf)
    df = pd.read_csv(f,sep='\t')
    # only take where Reviewed == reviewed
    r_df = df[df['Reviewed']=='reviewed'] 
    #  take columns ['Entry','Protein names','Gene names','Organism','EC number','Status']
    r_df = r_df[['Entry','Gene Names','EC number','Gene Ontology (GO)']]
    if os.path.exists(gbkf):
        print('gbkf exist:',gbkf)
        lt2ec,lt2oldlt, oldlt2lt = parse_genebank(gbkf)
    else:
        lt2ec = {}
        lt2oldlt = {}
        oldlt2lt = {}
    return r_df,df,lt2ec,lt2oldlt, oldlt2lt

def check_correct(d_list,lable):
    d =set()
    for i in d_list:
        for j in i.keys():
            d.add(j)
    for item in d:
        if item not in lable:
            return False
    return True


def get_labels_dis(file_name,pref,pred_names):
    gbkf = pref+'.gb'
    print('gbkf',gbkf)
    r_df,df,lt2ec,lt2oldlt, oldlt2lt = read_unif_genebk(file_name,gbkf)
    print('len of pred_names',len(pred_names))
    print('lt2ec',len(lt2ec.keys()),list(lt2ec.keys())[:10])
    ### remove lines with no EC number
    r_df = r_df.dropna(subset=['EC number'])
    print('Number of proteins with EC number:',len(r_df))
    all_label = set()
    label_dis={}
    if not lt2ec=={}:
        for i in range(0,len(pred_names)):
            protein = pred_names[i]
            try:
                oldname = lt2oldlt[protein]
            except:
                pass

            if len(oldname) ==1:
                # protein = oldname[0]
                unip  = r_df[r_df['Gene Names'].str.contains(oldname[0]+' ', na=False) | r_df['Gene Names'].str.endswith(oldname[0], na=False)]
            else:
                for old in oldname:
                    unip = r_df[r_df['Gene Names'].str.contains(old+' ', na=False) | r_df['Gene Names'].str.endswith(old, na=False)]
                    if not unip.empty:
                        # protein = old
                        break
            if unip.empty:
                continue
            else:
                ecnumber = unip['EC number'].values[0]
                ecs = ecnumber.split('; ')
                for ec in ecs:
                    label = ec.strip('EC:')
                    all_label.add(label)
                    if protein in label_dis:
                        label_dis[protein].update({label:1})
                    else:
                        label_dis[protein] = {label:1}
    else:
        for i in range(0,len(pred_names)):
            protein = pred_names[i]
            unip = r_df[r_df['Gene Names'].str.contains(protein+' ', na=False) | r_df['Gene Names'].str.endswith(protein, na=False)]
            if unip.empty:
                continue
            else:
                ecnumber = unip['EC number'].values[0]
                ecs = ecnumber.split('; ')
                for ec in ecs:
                    label = ec.strip('EC:')
                    all_label.add(label)
                    if protein in label_dis:
                        label_dis[protein].update({label:1})
                    else:
                        label_dis[protein] = {label:1}
    # vector_dis = label_dict_to_distribution(label_dis,all_label)
    return all_label,label_dis

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--file_name', type=str,required=False, default='/home/kexin/code/bigg/data/genomes/uniprot_269799.tsv', help='uniprot file')
parser.add_argument('--pref', type=str, default='/home/kexin/code/bigg/data/genomes/CP000148.1', help='prefix of genebank file')
parser.add_argument('--modelname', type=str, default='iAF987', help='model name')
parser.add_argument('--orif', type=str, default='../data/tmp/iAF987_allec_R1_T5_ORIpredscore_0.pkl', help='model name')
# parser.add_argument('--dig', type=int, default='/home/kexin/code/bigg/data/genomes/CP000148.1', help='prefix of genebank file')

args = parser.parse_args()
# file_name='/home/kexin/code/bigg/data/genomes/uniprot_269799.tsv'
# pref='/home/kexin/code/bigg/data/genomes/CP000148.1'
# modelname = 'iAF987'
# # cleandict = cleancsv2dic(pref+'_maxsep.csv')
# orif='/home/kexin/code/ibex_cobra/data/tmp/test_cp148_t5.0_ORIpredscore_0.pkl'

file_name=args.file_name
pref=args.pref
modelname = args.modelname
orif=args.orif
# modelnames=(iLJ478 iJN1463 iCN900 iAF1260 iAF987 iNJ661 iYO844)
# seqnames=(NC_000853.1_t1 NC_002947.4_t2 AM180355.1_t2 NC_000913.3_t4 CP000148.1_t4 NC_000962.3_t5 AL009126.3_t4)
# taxoids=(243274 160488 272563 83333 269799 83332 224308)
if modelname == 'iLJ478':
    pref = '/ibex/user/niuk0a/funcarve/cobra/uniprot/NC_000853.1'
    file_name = '/ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_243274.tsv'
elif modelname == 'iJN1463':
    pref = '/ibex/user/niuk0a/funcarve/cobra/uniprot/NC_002947.4'
    file_name = '/ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_160488.tsv'
elif modelname == 'iCN900':
    pref = '/ibex/user/niuk0a/funcarve/cobra/uniprot/AM180355.1'
    file_name = '/ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_272563.tsv'
elif modelname == 'iAF1260':
    pref = '/ibex/user/niuk0a/funcarve/cobra/uniprot/NC_000913.3'
    file_name = '/ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_83333.tsv'
elif modelname == 'iAF987':
    pref = '/ibex/user/niuk0a/funcarve/cobra/uniprot/CP000148.1'
    file_name = '/ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_269799.tsv'
elif modelname == 'iNJ661':
    pref = '/ibex/user/niuk0a/funcarve/cobra/uniprot/NC_000962.3'
    file_name = '/ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_83332.tsv'
elif modelname == 'iYO844':
    pref = '/ibex/user/niuk0a/funcarve/cobra/uniprot/AL009126.3'
    file_name = '/ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_224308.tsv'
modeltyp = orif.split('/')[-1].split('_')[1]
reward = orif.split('/')[-1].split('_')[2]
thr = orif.split('/')[-1].split('_')[3]
# file_name='/ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_269799.tsv'
# pref='/ibex/user/niuk0a/funcarve/cobra/uniprot/CP000148.1'
# modelname = 'iAF987'
# # cleandict = cleancsv2dic(pref+'_maxsep.csv')
# orif='../data/tmp/iAF987_allec_R1_T5_ORIpredscore_0.pkl'
# python check_turnover.py --file_name '' --pref '' --modelname iYO844 --orif ../data/tmp/iYO844_allec_R1_T5_ORIpredscore_0.pkl
with open(orif, 'rb') as file:
    oridata = pickle.load(file)
pred_names = list(oridata.keys())

tmpfile = orif.replace('ORIpredscore_0.pkl','differ_predscore')
dd={}
for i in range(1,4):
    f = tmpfile+'_'+str(i)+'.pkl'
    with open(f, 'rb') as file:
        data = pickle.load(file)
    dif = differ(oridata,data)
    for key in dif.keys():
        if key in dd:
            dd[key].append(dif[key])
        else:
            dd[key] = [dif[key]]
    break

       
all_label,label_dis = get_labels_dis(file_name,pref,pred_names)

# def check_correct_dig(d_list,lable,dig=3):
#     d =set()
#     for i in d_list:
#         for j in i.keys():
#             d.add(j)
#     lable = set(lable)
#     dlable = [x.split('.')[:dig] for x in lable]
#     # print(dlable)
#     for item in d:
#         ditem = item.split('.')[:dig]
#         if ditem not in dlable:
#             return False
#     return True

print('len of dd',len(dd))

cc = 0
total = 0
totalacc = 0
totalf1 = 0
totalprec = 0
totalrecall = 0
over=0 
overcc=0
print('dd.keys',list(dd.keys())[:10])
print('label_dis.keys',list(label_dis.keys())[:10])
for pr in dd.keys():
    if pr not in label_dis.keys():
        continue
    else:
        total+=1
        predecs =set()
        # avgval = sum(dd[pr].values())/len(dd[pr])
        # if avgval >= int(thr):
        #     over+=1
        avgval = 0
        for d in dd[pr]:
            predecs.update(d.keys())
            avgval += sum(d.values())/len(d)
            
        avgval = avgval/len(dd[pr])
        if avgval >= int(thr.replace('T','')):
            over+=1
        labelecs = set(label_dis[pr].keys())
        acc, recall, precison,f1 = get_metric(predecs,labelecs)
        # print(f'acc {acc}, f1 {f1}, precison {precison}, recall {recall}',flush=True)
        totalacc+=acc
        totalf1+=f1
        totalprec+=precison
        totalrecall+=recall
        if acc == 1:
            cc+=1
            overcc+=1
# print('name:',orif.split('/')[-1])
# # orif='../data/tmp/iAF987_allec_R1_T5_ORIpredscore_0.pkl'
# print(f'total {len(dd)} protein, in which {total} has label, {cc} full match, acc {totalacc/total:.3f}, f1 {totalf1/total:.3f}, precison {totalprec/total:.3f}, recall {totalrecall/total:.3f},over {over},overcc {overcc}',flush=True)
# exit()

with open('/ibex/user/niuk0a/funcarve/cobra/data/result/eval_turnover_7v4_I1.csv','a') as f:
    writer = csv.writer(f)
    # writer.writerow([modelname,modeltyp,reward,thr,len(dd),total,cc,total/len(dd),cc/total,totalacc/total,totalf1/total,totalprec/total,totalrecall/total])
    writer.writerow([modelname,modeltyp,reward,thr,len(dd),total,cc,total/len(dd),cc/total,totalacc/total,totalf1/total,totalprec/total,totalrecall/total,over,overcc])
print(f'total {len(dd)} protein, in which {total} has label, {cc} full match, acc {totalacc/total:.3f}, f1 {totalf1/total:.3f}, precison {totalprec/total:.3f}, recall {totalrecall/total:.3f}')
            
## iAF987_allec_R1_T5_differ_predscore_
# total 127 protein, in which 15 has label, 7 full match, acc 0.6533333333333334, f1 0.6981481481481481, precison 0.7166666666666667, recall 0.6866666666666668


# print(cc/total,total,len(dd),total/len(dd),cc)
# full match 0.5 12 102 0.11764705882352941 6
# dig3 0.75 12 102 0.11764705882352941 9
# dig2 0.8333333333333334 12 102 0.11764705882352941 10
    
