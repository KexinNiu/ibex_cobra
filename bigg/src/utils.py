import pandas as pd
import numpy as np
# from bigg.utils import *
import pickle
# import cobra
# import os
# print(os.getcwd())

def read_rhea2go(f:str):
    rhea2go = {}
    with open(f, 'r') as f:
        for line in f:
            if line.startswith('!'):
                continue    
            line = line.strip('\n').split(' ')
            rhea2go[line[0]] = line[-1]
    return rhea2go

def get_rhea_go(rhea2go:dict, rhea_ids:list):
    go = set()
    for i in rhea_ids:
        if i in rhea2go:
            go.add(rhea2go[i])
    return go

def read_ec2go(f:str):
    ec2go = {}
    with open(f, 'r') as f:
        for line in f:
            if line.startswith('!'):
                continue    
            line = line.strip('\n').split(' ')
            ec2go[line[0]] = line[-1]
    return ec2go

def get_ec_go(ec2go:dict, ec:list):
    go = set()
    for i in ec:
        if i in ec2go:
            go.add(ec2go[i])
    return go

def read_keggr2go(f:str):
    keggr2go = {}
    with open(f, 'r') as f:
        for line in f:
            if line.startswith('!'):
                continue    
            line = line.strip('\n').split(' ')
            keggr2go[line[0]] = line[-1]
    return keggr2go

def get_keggr_go(keggr2go:dict, keggr:list):
    go = set()
    for i in keggr:
        if i in keggr2go:
            go.add(keggr2go[i])
    return go

def database_links_reformat(database_links:list):
    links = []
    # print('database_links',database_links)
    if type(database_links) == str:
        database_links = [database_links]
    if type(database_links) != list:
        # print('database_links',database_links)
        return []
    for i in database_links:
        if ';' in i:
            i = i.split('; ')
            for j in i:
                links.append(j)
        else:
            links.append(i)
    return links 

def links_to_id(links:list):
    if links == []:
        return [], [], [], [], [], [], []
    ids = []
    rhea_ids = []
    mnxs = []
    seeds = []
    biocyc = []
    ec = []
    keggr = []
    for i in links:
        if i.startswith('RHEA'):
            id = i.split('/')[-1]
            rhea_ids.append(f'RHEA:{id}')
            ids.append(f'RHEA:{id}')
        elif i.startswith('MetaNetX'):
            id = i.split('/')[-1]
            mnxs.append(f'MNX:{id}')
            ids.append(f'MNX:{id}')
        elif i.startswith('SEED'):
            id = i.split(':')[-1]
            seeds.append(f'SEED:{id}')
            ids.append(f'SEED:{id}')
        elif i.startswith('BioCyc'):
            id = i.split('/')[-1]
            biocyc.append(f'BioCyc:{id}')
            ids.append(f'BioCyc:{id}')
        elif i.startswith('EC Number:'):
            id = i.split('/')[-1]
            ec.append(f'EC:{id}')
            ids.append(f'EC:{id}')
        elif i.startswith('KEGG Reaction'):
            id = i.split('/')[-1]
            keggr.append(f'KEGG_REACTION:{id}')
            ids.append(f'KEGG_REACTION:{id}')
        # RHEA:001
        # MNX:001
        # SEED:001
    return ids, rhea_ids, mnxs, seeds, biocyc, ec, keggr




def read_fasta(fasta_file):
    seq=''
    names =[]
    seqs = []
    with open(fasta_file, 'r') as inFile:
        for line in inFile:
            if line.startswith('>'):
                try:
                    name = line.strip('\n').split('[gene=')[1].split(']')[0]
                except IndexError:
                    name = line.strip('\n').split('[locus_tag=')[1].split(']')[0]   
                names.append(name)
                if seq == '':
                    continue
                else:
                    seqs.append(seq)
                    seq = ''    
            else:
                seq = seq + line.strip('\n')
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
    
def read_clean_withscore(input_file):
    threrhold = 0.8
    print('threrhold-->',threrhold)
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
                    if dis >= threrhold:
                        try:
                            pr2ec[pr].append(ecid)
                            # predscore[pr].update({ecid:dis})
                        except KeyError:
                            pr2ec[pr] = [ecid]   
                            # predscore[pr] = {ecid:dis} 
    print('pr2ec-protein number->',len(list(pr2ec.keys())))   
    return pr2ec,predscore

def clean2biggr(universal,predscore):
    print('universal',universal)
    universal_scoredict={}
    u_rx = [r.id for r in universal.reactions]
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

seednum = {
    0: '555',
    1: '666',
    2: '777',
    3: '888',
    4: '999'
}

media = 'default'
modelseed_to_bigg = {
    'cpd00001_e': 'h2o_e', 'cpd00035_e': 'o2_e', 'cpd00041_e': 'co2_e', 'cpd00023_e': 'glc__D_e',
    'cpd00119_e': 'nh4_e', 'cpd00107_e': 'pi_e', 'cpd00060_e': 'so4_e', 'cpd00161_e': 'k_e',
    'cpd00069_e': 'fe2_e', 'cpd00084_e': 'h_e', 'cpd00033_e': 'ac_e', 'cpd00322_e': 'cl_e',
    'cpd00066_e': 'mg2_e', 'cpd00054_e': 'na1_e', 'cpd00065_e': 'ca2_e', 'cpd00156_e': 'cu2_e',
    'cpd00220_e': 'mn2_e', 'cpd00644_e': 'zn2_e', 'cpd00393_e': 'cobalt2_e', 'cpd00133_e': 'ni2_e',
    'cpd00263_e': 'mobd_e', 'cpd00104_e': 'trp__L_e', 'cpd00149_e': 'his__L_e', 'cpd00971_e': 'gly_e',
    'cpd00099_e': 'ala__L_e', 'cpd00205_e': 'ser__L_e', 'cpd00009_e': 'nad_e', 'cpd00063_e': 'asp__L_e',
    'cpd00254_e': 'glu__L_e', 'cpd10515_e': 'phe__L_e', 'cpd00030_e': 'arg__L_e', 'cpd00242_e': 'leu__L_e',
    'cpd00226_e': 'ile__L_e', 'cpd01242_e': 'val__L_e', 'cpd00307_e': 'thr__L_e', 'cpd00092_e': 'lys__L_e',
    'cpd00117_e': 'met__L_e', 'cpd00067_e': 'pro__L_e', 'cpd00567_e': 'cys__L_e', 'cpd00132_e': 'asn__L_e',
    'cpd00210_e': 'gln__L_e', 'cpd00320_e': 'tyr__L_e', 'cpd03279_e': 'orn_e', 'cpd00246_e': 'tryptamine_e',
    'cpd00311_e': 'pyr_e', 'cpd00367_e': 'glycine_e', 'cpd00277_e': 'phenylalanine_e', 'cpd00182_e': 'oxaloacetate_e',
    'cpd00654_e': 'oxalate_e', 'cpd00412_e': 'fumarate_e', 'cpd00438_e': 'succinate_e', 'cpd00274_e': 'malate_e',
    'cpd00186_e': 'citrate_e', 'cpd00637_e': 'lactate_e', 'cpd00105_e': 'aspartate_e', 'cpd00305_e': 'alpha_KG_e',
    'cpd00309_e': 'pyruvate_e', 'cpd00098_e': 'glutamate_e', 'cpd00207_e': 'glutamine_e', 'cpd00082_e': 'formate_e',
    'cpd00129_e': 'succ_e'
}
# Load biggr2ec dictionary from the file
with open('../data/biggr2ec.pkl', 'rb') as f:
    biggr2ec = pickle.load(f)
with open('../data/biggec2r.pkl', 'rb') as f:
    biggec2r = pickle.load(f)
    # cobra/bigg/data/biggec2r.pkl
print("biggr2ec dictionary loaded successfully.")
# Your media setup code
# if media == 'rich':
#     media = ['cpd00001_e','cpd00035_e','cpd00041_e','cpd00023_e','cpd00119_e','cpd00107_e','cpd00060_e',
#              'cpd00161_e','cpd00069_e','cpd00084_e','cpd00033_e', 'cpd00322_e','cpd00066_e','cpd00054_e',
#              'cpd00065_e','cpd00156_e','cpd00220_e','cpd00644_e','cpd00393_e','cpd00133_e','cpd00263_e',
#              'cpd00104_e','cpd00149_e','cpd00971_e','cpd00099_e','cpd00205_e','cpd00009_e','cpd00063_e',
#              'cpd00254_e','cpd10515_e','cpd00030_e','cpd00242_e','cpd00226_e','cpd01242_e','cpd00307_e',
#              'cpd00092_e','cpd00117_e','cpd00067_e','cpd00567_e','cpd00132_e','cpd00210_e','cpd00320_e',
#              'cpd03279_e','cpd00246_e','cpd00311_e','cpd00367_e','cpd00277_e','cpd00182_e','cpd00654_e',
#              'cpd00412_e','cpd00438_e','cpd00274_e','cpd00186_e','cpd00637_e','cpd00105_e','cpd00305_e',
#              'cpd00309_e','cpd00098_e','cpd00207_e','cpd00082_e','cpd00129_e']
# elif media == 'minimal':
#     media = ['cpd00001_e','cpd00065_e','cpd00060_e','cpd00322_e','cpd00129_e','cpd00156_e','cpd00107_e',
#              'cpd00084_e','cpd00149_e','cpd00099_e','cpd10515_e','cpd00030_e','cpd00254_e','cpd00063_e',
#              'cpd00205_e','cpd00009_e','cpd00971_e','cpd00242_e','cpd00104_e','cpd00644_e','cpd00263_e',
#              'cpd00082_e']
# elif media == 'default':
#     media = ['cpd00035_e','cpd00051_e','cpd00132_e','cpd00041_e','cpd00084_e','cpd00053_e','cpd00023_e',
#              'cpd00033_e','cpd00119_e','cpd00322_e','cpd00107_e','cpd00039_e','cpd00060_e','cpd00066_e',
#              'cpd00129_e','cpd00054_e','cpd00161_e','cpd00065_e','cpd00069_e','cpd00156_e','cpd00027_e',
#              'cpd00149_e','cpd00030_e','cpd00254_e','cpd00971_e','cpd00063_e','cpd10515_e','cpd00205_e',
#              'cpd00099_e']
# else:
#     media = media
    
# media = [modelseed_to_bigg[cpd] if cpd in modelseed_to_bigg else cpd for cpd in media]