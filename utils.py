
import pandas as pd
import numpy as np
from bigg.utils import *
import pickle
import cobra
# modelseed_to_bigg = {
#     'cpd00001_e': 'h2o_e',
#     'cpd00035_e': 'o2_e',
#     'cpd00041_e': 'co2_e',
#     'cpd00023_e': 'glc__D_e',
#     'cpd00119_e': 'nh4_e',
#     'cpd00107_e': 'pi_e',
#     'cpd00060_e': 'so4_e',
#     'cpd00161_e': 'k_e',
#     'cpd00069_e': 'fe2_e',
#     'cpd00084_e': 'h_e',
#     'cpd00033_e': 'ac_e',
#     'cpd00322_e': 'cl_e',
#     'cpd00066_e': 'mg2_e',
#     'cpd00054_e': 'na1_e',
#     'cpd00065_e': 'ca2_e',
#     'cpd00156_e': 'cu2_e',
#     'cpd00220_e': 'mn2_e',
#     'cpd00644_e': 'zn2_e',
#     'cpd00393_e': 'cobalt2_e',
#     'cpd00133_e': 'ni2_e',
#     'cpd00263_e': 'mobd_e',
#     'cpd00104_e': 'trp__L_e',
#     'cpd00149_e': 'his__L_e',
#     'cpd00971_e': 'gly_e',
#     'cpd00099_e': 'ala__L_e',
#     'cpd00205_e': 'ser__L_e',
#     'cpd00009_e': 'nad_e',
#     'cpd00063_e': 'asp__L_e',
#     'cpd00254_e': 'glu__L_e',
#     'cpd10515_e': 'phe__L_e',
#     'cpd00030_e': 'arg__L_e',
#     'cpd00242_e': 'leu__L_e',
#     'cpd00226_e': 'ile__L_e',
#     'cpd01242_e': 'val__L_e',
#     'cpd00307_e': 'thr__L_e',
#     'cpd00092_e': 'lys__L_e',
#     'cpd00117_e': 'met__L_e',
#     'cpd00067_e': 'pro__L_e',
#     'cpd00567_e': 'cys__L_e',
#     'cpd00132_e': 'asn__L_e',
#     'cpd00210_e': 'gln__L_e',
#     'cpd00320_e': 'tyr__L_e',
#     'cpd03279_e': 'orn_e',
#     'cpd00246_e': 'tryptamine_e',
#     'cpd00311_e': 'pyr_e',
#     'cpd00051_e':
#     'cpd00367_e': 'glycine_e',
#     'cpd00277_e': 'phenylalanine_e',
#     'cpd00182_e': 'oxaloacetate_e',
#     'cpd00654_e': 'oxalate_e',
#     'cpd00412_e': 'fumarate_e',
#     'cpd00438_e': 'succinate_e',
#     'cpd00274_e': 'malate_e',
#     'cpd00186_e': 'citrate_e',
#     'cpd00637_e': 'lactate_e',
#     'cpd00105_e': 'aspartate_e',
#     'cpd00305_e': 'alpha_KG_e',
#     'cpd00309_e': 'pyruvate_e',
#     'cpd00098_e': 'glutamate_e',
#     'cpd00207_e': 'glutamine_e',
#     'cpd00082_e': 'formate_e',
#     'cpd00129_e': 'succ_e'
# }

def read_fasta(fasta_file):
    seq=''
    names =[]
    seqs = []
    with open(fasta_file, 'r') as inFile:
        for line in inFile:
            if line.startswith('>'):
                name = line.strip('\n').split('>')[1]
                # try:
                #     name = line.strip('\n').split('[gene=')[1].split(']')[0]
                # except IndexError:
                #     name = line.strip('\n').split('[locus_tag=')[1].split(']')[0]   
                names.append(name)
                if seq == '':
                    continue
                else:
                    seqs.append(seq)
                    seq = ''    
            else:
                seq = seq + line.strip('\n')
        seqs.append(seq)
    return names,seqs
    
            
def read_clean(input_file):
    threrhold = 0.8
    print('threrhold-->',threrhold)
    pr2ec = {}
    with open(input_file, 'r') as inFile:
        for line in inFile:
            line = line.split()
            pr = line[0]
            items = line[-1].split(',')
            for item in items:
                if item.startswith('EC:'):
                    ec,dis = item.split('/')
                    ecid = ec.split(':')[-1]
                    dis= float(dis)
                    if dis >= threrhold:
                        try:
                            pr2ec[pr].append(ecid)
                        except KeyError:
                            pr2ec[pr] = [ecid]     
    print('pr2ec-protein number->',len(list(pr2ec.keys())))   
    
def read_clean_withscore(input_file,threshold=0.8):
    print('threrhold-->',threshold)
    pr2ec = {}
    predscore = {}
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
                    # ecid = ec.split(':')[-1]
                    ecid = ec
                    dis = float(dis)
                    if dis >= 0.0001:
                        try:
                            predscore[pr].update({ecid:dis})
                        except:
                            predscore[pr] = {ecid:dis}
                    if dis >= threshold:
                        try:
                            pr2ec[pr].append(ecid)
                            # predscore[pr].update({ecid:dis})
                        except KeyError:
                            pr2ec[pr] = [ecid]   
                            # predscore[pr] = {ecid:dis} 
    print('pr2ec-protein number->',len(list(pr2ec.keys())))   
    return pr2ec,predscore


# # generate biggr2ec
# reactionf = '/ibex/user/niuk0a/funcarve/cobra/bigg/bigg_models_reactions.txt'
# reactions = pd.read_csv(reactionf, sep='\t')
# ecs = []

# for link in list(reactions.database_links.values):
#     link = database_links_reformat(link)
#     ids, rhea_ids, mnxs, seeds, biocyc, ec, keggr = links_to_id(link)
#     ecs.append(ec)

# biggr2ec = {}
# biggec2r = {}
# for index, row in reactions.iterrows():
#     id = row['bigg_id']
#     oldid = row['old_bigg_ids']
    
#     ec = ecs[index]
#     for eitem in ec : 
#         try:
#             biggec2r[eitem].append(id)
#             if '; ' in oldid:
#                 oids = oldid.split('; ')
#                 for oid in oids:
#                     biggec2r[eitem].append(oid) 
#             else:
#                 biggec2r[eitem].append(oldid)
#         except:
#             biggec2r[eitem] = [id]
#             if '; ' in oldid:
#                 oids = oldid.split('; ')
#                 biggec2r[eitem] = oids
#             else:
#                 biggec2r[eitem] = [oldid]
#     try:
#         biggr2ec[id] += ec
#         if '; ' in oldid:
#             oids = oldid.split('; ')
#             for oid in oids:
#                 biggr2ec[oid] += ec
#         else:
#             biggr2ec[oldid] += ec
#     except:
#         biggr2ec[id] = ec
#         if '; ' in oldid:
#             oids = oldid.split('; ')
#             for oid in oids:
#                 biggr2ec[oid] = ec
#         else:
#             biggr2ec[oldid] = ec
       
# #save biggr2ec


# # Save biggr2ec dictionary to a file
# with open('biggr2ec.pkl', 'wb') as f:
#     pickle.dump(biggr2ec, f)
# with open('biggec2r.pkl', 'wb') as f:
#     pickle.dump(biggec2r, f)
# print("biggr2ec dictionary saved successfully.")


# Load biggr2ec dictionary from the file
with open('biggr2ec.pkl', 'rb') as f:
    biggr2ec = pickle.load(f)
with open('biggec2r.pkl', 'rb') as f:
    biggec2r = pickle.load(f)
print("biggr2ec dictionary loaded successfully.")



def clean2biggr(predscore):
    # print('universal',universal)
    universal_scoredict={}
    # u_rx = [r.id for r in universal.reactions]
    for pr, ec_score in predscore.items():
        for ec , score in ec_score.items():
            if ec in biggec2r.keys():
                rids = biggec2r[ec]
                # if len(rids) >1:
                #     print('rids:',rids)
                for rid in rids:
                    # print('rid:',rid)
                    try:
                        universal_scoredict[rid].append(score)
                    except:
                        universal_scoredict[rid] = score

    return universal_scoredict

def clean_to_rxns(pr2ec,r2ecf,pr2gene,gene_modelseed,organism):
    # loose =True
    loose = False
    ecf = pd.read_csv(r2ecf, sep='\t')
    print('done read ecf',flush=True)
    if loose:
        allshort = set()
        for ec in ecf['External ID'].values:
            if '.-' in ec and ec.count('-') <= 2:
                allshort.add(ec)

    def check(ecid,allshort):
        ecid1 = ecid.split('.')[:-1]
        ecid1 = '.'.join(ecid1)+'.-'
        ecid2 = ecid.split('.')[:-2]
        ecid2 = '.'.join(ecid2)+'.-.-'
        if ecid1 in allshort:
            return True,ecid1
        elif ecid2 in allshort:
            return True,ecid2
        else:
            return False,0
        
    print('if loose is true, we will add the loose reactions,loose=',loose,flush=True)
    
    # if org != 'default':
    #     new_hits = _get_org_rxns(gene_modelseed, organism)
    #     gene_count = len(new_hits)
    #     print('Added', gene_count, 'genes from', organism,flush=True)
    rxn_p = {}
    print('using origianl enzyme to reaction mapping')
    for pr in pr2ec.keys():
        ecs = pr2ec[pr]
        # gene = pr 
        # try:
        #     gene = pr2gene[pr]
        # except KeyError:
        #     continue
        # try:
        #     rxns = gene_modelseed[gene]
        # except KeyError:
        #     continue
        for ecid in ecs:
            rs = ecf.loc[(ecf['External ID'] == ecid), 'ModelSEED ID'].values.tolist()
            if loose:
                f,cutid =check(ecid,allshort)
                if f:
                    loose_rxns = ecf.loc[(ecf['External ID'] == cutid), 'ModelSEED ID'].values.tolist()
                    rs = rs + loose_rxns
                    rs = list(set(rs))
            else:
                rs = list(set(rs))
            for r in rs:
                r = r + '_c'
                try:
                    rxn_p[r].append(pr)
                except KeyError:
                    rxn_p[r] = [pr]
                
    for r in rxn_p.keys():
        rxn_p[r] = list(set(rxn_p[r]))
    print('rxn_p rxn number:',len(list(rxn_p.keys())),flush=True)
    return rxn_p

