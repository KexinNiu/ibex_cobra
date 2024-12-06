import pandas as pd
import csv
import numpy as np
import os
from Bio import SeqIO
import argparse


def softmax(vector):
    # Add a small constant to avoid division by zero and ensure numerical stability
    exp_vector = np.exp(vector - np.max(vector))  # Subtract max for numerical stability
    sum_exp = np.sum(exp_vector)
    
    # Check if the sum is zero (this would mean the vector has all zeros)
    if sum_exp == 0:
        return np.zeros_like(vector)  # Return a zero vector or use uniform distribution
    
    # Normalize the exponentiated vector
    return exp_vector / sum_exp


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
       
def read_unif(f):
    df = pd.read_csv(f,sep='\t')
    # only take where Reviewed == reviewed
    r_df = df[df['Reviewed']=='reviewed'] 
    #  take columns ['Entry','Protein names','Gene names','Organism','EC number','Status']
    r_df = r_df[['Entry','Gene Names','EC number','Gene Ontology (GO)']]
    return r_df,df


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

import numpy as np

def kl_divergence(P, Q, epsilon=1e-10):
    # Normalize P and Q to make them valid probability distributions
    # P = P / np.sum(P)
    # Q = Q / np.sum(Q)
    
    # # Add a small value to avoid log(0) issues
    # P = np.where(P == 0, epsilon, P)
    # Q = np.where(Q == 0, epsilon, Q)
    P = np.clip(P, epsilon, 1)  # Ensure P does not have zero values
    Q = np.clip(Q, epsilon, 1) 
    # Calculate KL divergence
    return np.sum(P * np.log(P / Q))

# # Example vectors
# P = np.array([0.2, 0.5, 0.3])
# Q = np.array([0.1, 0.4, 0.5])

# kl_div = kl_divergence(P, Q)
# print("KL Divergence:", kl_div)
def cleancsv2dic(f):
    with open(f, 'r') as file:
        reader = csv.reader(file)
        d = {}
        for row in reader:
            name = row[0]
            ecs = row[1:]
            for ec in ecs:
                label,score = ec.split('/')
                label = label.strip('EC:')
                score = float(score)
                if name in d:
                    d[name].update({label:score})
                else:
                    d[name] = {label:score}
    return d

import itertools


def label_dict_to_distribution(protein_labels, all_label):
    # Create a list of all possible labels and sort to maintain order consistency
    all_label = sorted(all_label)
    label_indices = {label: idx for idx, label in enumerate(all_label)}
    
    # Initialize the result dictionary
    result = {}
    result_softmax = {}

    # Iterate over each protein and its associated labels and scores
    for protein, labels_scores in protein_labels.items():
        labels = list(labels_scores.keys())  # Extract the labels for the protein
        
        # Generate all non-empty combinations of labels
        label_combinations = []
        for r in range(1, len(labels) + 1):
            label_combinations.extend(itertools.combinations(labels, r))
        
        # Convert each combination to a one-hot vector based on `all_label`
        protein_vectors = []
        for combination in label_combinations:
            vector = np.zeros(len(all_label))
            for label in combination:
                vector[label_indices[label]] = labels_scores[label]  # Set the score for this label
            protein_vectors.append(vector)
        
        # softmax the vectors
        prv_softmax = [softmax(vector) for vector in protein_vectors]
        # Store the list of vectors for this protein in the result dictionary
        result[protein] = protein_vectors
        result_softmax[protein] = prv_softmax
    
    return result


def label_dict_to_one_hot(protein_labels, all_label):
    # Create a sorted list of all possible labels and a mapping to indices
    all_label = sorted(all_label)
    label_indices = {label: idx for idx, label in enumerate(all_label)}
    
    result = {}
    
    # Process each protein and its associated labels
    for protein, labels_scores in protein_labels.items():
        # Initialize the one-hot vector with zeros
        one_hot_vector = np.zeros(len(all_label))
        
        # For each label in the protein_labels, set the corresponding index in the one-hot vector
        for label, score in labels_scores.items():
            if label in label_indices:  # Ensure the label is in all_label
                one_hot_vector[label_indices[label]] = score
        
        # Handle the case where the one-hot vector is all zeros (no labels for the protein)
        if np.sum(one_hot_vector) == 0:
            # If no labels, set softmax to uniform distribution or zero
            softmax_vector = np.zeros(len(all_label))
        else:
            # Apply softmax to the one-hot vector
            softmax_vector = softmax(one_hot_vector)
        
        # Store the softmax-processed vector in the result dictionary
        result[protein] = softmax_vector
    
    return result

def get_labels_dis(file_name,pref,pred_names):
    # pref = pref.replace('/ibex/user/niuk0a/CLEAN/app/results/inputs/','/ibex/user/niuk0a/funcarve/cobra/uniprot/')
    # /ibex/user/niuk0a/funcarve/cobra/uniprot/NC_000853.1.gb
    gbkf = pref+'.gb'
    # print('gene bank file:',gbkf)
    r_df,df,lt2ec,lt2oldlt, oldlt2lt = read_unif_genebk(file_name,gbkf)
    ### remove lines with no EC number
    r_df = r_df.dropna(subset=['EC number'])
    # print('r_df:',r_df.head())
    print('Number of proteins with EC number:',len(r_df))
    all_label = set()
    label_dis={}
    if not lt2oldlt=={}:
        print('Using locus tag to EC number mapping')
        for i in range(0,len(pred_names)):
            protein = pred_names[i]
            try:
                oldname = lt2oldlt[protein]
            except:
                pass

            for old in oldname:
                # print('old:',old)
                unip = r_df[r_df['Gene Names'].str.contains(old+' ', na=False) | r_df['Gene Names'].str.endswith(old, na=False)]
                # protein = old
                if not unip.empty:
                    ecnumber = unip['EC number'].values[0]
                    ecs = ecnumber.split('; ')
                    for ec in ecs:
                        label = ec.strip('EC:')
                        all_label.add(label)
                        if old in label_dis:
                            label_dis[protein].update({label:1})
                        else:
                            label_dis[protein] = {label:1}
    else:
        print('not using lt2ec')
        for i in range(0,len(pred_names)):
            protein = pred_names[i]
            unip = r_df[r_df['Gene Names'].str.contains(protein+' ', na=False) | r_df['Gene Names'].str.endswith(protein, na=False)]
            if unip.empty:
                # print('empty:',protein)
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
    # print('label dis:',len(label_dis),list(label_dis.keys())[:20])
    vector_dis = label_dict_to_distribution(label_dis,all_label)
    return all_label,label_dis,vector_dis

def genome_kl(predvector,truevector):
    kl = 0
    for protein in truevector.keys():
        # if protein not in predvector.keys():
        #     # print('protein not in predvector:',protein)
        #     continue
        # print(protein)
        # print(len(predvector[protein]))
        # print(len(truevector[protein][0]))
        ## add the minimum kl divergence between the predicted and one of the true distribution
        score = [kl_divergence(predvector[protein],truevector[protein][i]) for i in range(len(truevector[protein]))]
        kl += min(score)
        # break
    return kl,kl/len(truevector.keys())

parser = argparse.ArgumentParser()
parser.add_argument('--file_name', type=str, help='file name')
parser.add_argument('--pref', type=str, help='prefix')
parser.add_argument('--modelname', type=str, help='model name')
parser.add_argument('--index', type=int, help='index')
args = parser.parse_args()
file_name = args.file_name
modelname = args.modelname
index = args.index
pref = args.pref



# file_name='/home/kexin/code/bigg/data/genomes/uniprot_269799.tsv'
# pref='/home/kexin/code/bigg/data/genomes/CP000148.1'
# modelname = 'iAF987'
# python 00_eval.py --file_name /ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_243274.tsv --pref /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1  --modelname iLJ478 --index 1
# file_name = '/ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_243274.tsv'
# pref = '/ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1'
# modelname = 'iLJ478'
# index = 1




print('clean file:',pref+'_maxsep.csv')
cleandict = cleancsv2dic(pref+'_maxsep.csv')
pred_names = list(cleandict.keys())
# print(len(pred_names))
# print(len(all_label))
# print(len(label_dis))
# print('vector_dis',len(vector_dis),list(vector_dis.keys())[:20])
all_label,label_dis,vector_dis = get_labels_dis(file_name,pref,pred_names)

cleanvector_dis = label_dict_to_one_hot(cleandict,all_label)
# # print('cleanvector_dis',len(cleanvector_dis),list(cleanvector_dis.keys())[:20])
cleankl,cleanavgkl = genome_kl(cleanvector_dis,vector_dis)
tmp = 'CLEAN KL|'+str(cleankl)+'|'+str(cleanavgkl)


result = [tmp]
# result = []
print('clean t files:',pref+f'_t{index}_maxsep.csv')
cleandict = cleancsv2dic(pref+f'_t{index}_maxsep.csv')
pred_names = list(cleandict.keys())
print(len(pred_names))
all_label,label_dis,vector_dis = get_labels_dis(file_name,pref,pred_names)
cleanvector_dis = label_dict_to_one_hot(cleandict,all_label)
cleankl,cleanavgkl = genome_kl(cleanvector_dis,vector_dis)
tmp = 'CLEAN KLt4|'+str(cleankl)+'|'+str(cleanavgkl)
result.append(tmp)
# newprediscoref = '/home/kexin/code/bigg/data/genomes/results/' + str(modelname) +f'_Rt8_t0_newpredscore_'+str(i)+'.pkl'
# newprediscoref = '/home/kexin/code/bigg/data/genomes/results/iAF987_Rt8_t0_newpredscore_9.pkl'
for i in range(1,10):
# for i in [1,5,9]:
    # newprediscoref = '/home/kexin/code/bigg/data/genomes/results/' + str(modelname) +f'_R_t{thr}_newpredscore_'+str(i)+'.pkl'
    print(i)
    # /ibex/user/niuk0a/funcarve/cobra/iAF987_Rt8_t0_newpredscore_1.pkl
    # /ibex/user/niuk0a/funcarve/cobra/iAF987_Rt8_t0_newpredscore_1.pkl
    if modelname == 'iAF987':
        newprediscoref = '/ibex/user/niuk0a/funcarve/cobra/' + str(modelname) +f'_Rt8_t0_newpredscore_'+str(i)+'.pkl'
    else:
        newprediscoref = '/ibex/user/niuk0a/funcarve/cobra/' + str(modelname) +f'_Rt{index}_t8_newpredscore_'+str(i)+'.pkl'
    if not os.path.exists(newprediscoref):
        break
    # newprediscoref = f'/home/kexin/code/bigg/data/genomes/results/iAF987_R_t5_updateprs_{i}.pkl'
    interdict = pd.read_pickle(newprediscoref)
    # print(interdict,len(interdict.keys()))
    intervector_dis = label_dict_to_one_hot(interdict,all_label)
    interkl,interavgkl = genome_kl(intervector_dis,vector_dis)
    tmp = f'INTER {i}|'+str(interkl)+'|'+str(interavgkl)
    result.append(tmp)
    
    
print(modelname)
for r in result:
    print(r)




