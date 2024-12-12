import pickle
import pandas as pd
import csv
import numpy as np
import os
from Bio import SeqIO
import argparse




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
    df = pd.read_csv(f,sep='\t')
    # only take where Reviewed == reviewed
    r_df = df[df['Reviewed']=='reviewed'] 
    #  take columns ['Entry','Protein names','Gene names','Organism','EC number','Status']
    r_df = r_df[['Entry','Gene Names','EC number','Gene Ontology (GO)']]
    if os.path.exists(gbkf):
        lt2ec,lt2oldlt, oldlt2lt = parse_genebank(gbkf)
    else:
        lt2ec = {}
        lt2oldlt = {}
        oldlt2lt = {}
    return r_df,df,lt2ec,lt2oldlt, oldlt2lt


def get_labels_dis(file_name,pref,pred_names):
    gbkf = pref+'.gb'
    r_df,df,lt2ec,lt2oldlt, oldlt2lt = read_unif_genebk(file_name,gbkf)
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
                protein = oldname[0]
                unip  = r_df[r_df['Gene Names'].str.contains(oldname[0]+' ', na=False) | r_df['Gene Names'].str.endswith(oldname[0], na=False)]
            else:
                for old in oldname:
                    unip = r_df[r_df['Gene Names'].str.contains(old+' ', na=False) | r_df['Gene Names'].str.endswith(old, na=False)]
                    if not unip.empty:
                        protein = old
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
parser.add_argument('--file_name', type=str, default='/home/kexin/code/bigg/data/genomes/uniprot_269799.tsv', help='uniprot file')
parser.add_argument('--pref', type=str, default='/home/kexin/code/bigg/data/genomes/CP000148.1', help='prefix of genebank file')
parser.add_argument('--modelname', type=str, default='iAF987', help='model name')
parser.add_argument('--dig', type=int, default='/home/kexin/code/bigg/data/genomes/CP000148.1', help='prefix of genebank file')

args = parser.parse_args()
file_name='/home/kexin/code/bigg/data/genomes/uniprot_269799.tsv'
pref='/home/kexin/code/bigg/data/genomes/CP000148.1'
modelname = 'iAF987'
# cleandict = cleancsv2dic(pref+'_maxsep.csv')
orif='/home/kexin/code/ibex_cobra/data/tmp/test_cp148_t5.0_ORIpredscore_0.pkl'
with open(orif, 'rb') as file:
    oridata = pickle.load(file)
pred_names = list(oridata.keys())

all_label,label_dis = get_labels_dis(file_name,pref,pred_names)

dd={}
for i in range(1,4):
    f =f'/home/kexin/code/ibex_cobra/data/tmp/test_cp148_t5_differ_predscore_{i}.pkl'
    # f = f'test_cp148_t5_updateprs_{i}.pkl'
    with open(f, 'rb') as file:
        data = pickle.load(file)
    for key in data.keys():
        if key not in dd:
            dd[key] = []
        dd[key].append(data[key])

def check_correct(d_list,lable):
    d =set()
    for i in d_list:
        for j in i.keys():
            d.add(j)
    for item in d:
        if item not in lable:
            return False
    return True
def check_correct_dig(d_list,lable,dig=3):
    d =set()
    for i in d_list:
        for j in i.keys():
            d.add(j)
    lable = set(lable)
    dlable = [x.split('.')[:dig] for x in lable]
    # print(dlable)
    for item in d:
        ditem = item.split('.')[:dig]
        if ditem not in dlable:
            return False
    return True


cc = 0
total = 0
for key in dd.keys():
    print(key,end='|')
    for i in range(len(dd[key])):
        print(dd[key][i],end='|')
    try:
        print(label_dis[key],end='++')
        total+=1
        if check_correct_dig(list(dd[key]),label_dis[key],2):
        # if list(dd[key][0].keys()) == list(label_dis[key].keys()) or list(dd[key][0].keys())[0] in list(label_dis[key].keys()):
            cc+=1
            print('|||equal')
        else:
            print(dd[key][-1],'|||not equal')
    except:
        print('no label')
        pass
   
print(cc/total,total,len(dd),total/len(dd),cc)
# full match 0.5 12 102 0.11764705882352941 6
# dig3 0.75 12 102 0.11764705882352941 9
# dig2 0.8333333333333334 12 102 0.11764705882352941 10
    
