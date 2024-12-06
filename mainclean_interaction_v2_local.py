#!/usr/bin/env

### using all possible ec from clean 
'''Reconstructor 

Reconstructor is an automatic genome scale metabolic network reconstruction tool that is user-friendly, COBRApy compatible, and uses a pFBA-based
gap-filling technique. 

Inputs
---------
Type 1: Annotated protein fasta file
Type 2: BLASTp output
Type 3: SBML Model

Output
---------
Well-annotated SBML model that uses ModelSEED namespace and is directly compatible with COBRApy withotut the need for additional compatibility modules

Example of how to run reconstructor
-----------------------------------
Type 1 input:  python -m reconstructor --input Osplanchnicus.aa.fasta --type 1 --gram negative --other_args <args>
Type 2 input: python -m reconstructor --input Osplanchnicus.hits.out --type 2 --gram negative --other_args <args>
Type 3 input: python -m reconstructor --input Osplanchnicus.sbml --type 3 --other_args <args>

Options for Running Reconstructor 
---------------------------------
--input <input file, Required>
--type <input file type, .fasta = 1, diamond blastp output = 2, .sbml = 3, Required, Default = 1> 
--gram <Type of Gram classificiation (positive or negative), default = positive>
--media <List of metabolites composing the media condition. Not required.>
--tasks <List of metabolic tasks. Not required>
--org <KEGG organism code. Not required>
--min_frac <Minimum objective fraction required during gapfilling, default = 0.01>
--max_frac <Maximum objective fraction allowed during gapfilling, default = 0.5>
--out <Name of output GENRE file, default = default>
--name <ID of output GENRE, default = default>
--cpu <Number of processors to use, default = 1>
--test <run installation tests, default = no>
'''

# Dependencies
import sys
import time
import numpy as np
import wget
import shutil
import os
import cobra
import pickle
import argparse
import warnings
import symengine
import subprocess
from random import shuffle
from multiprocessing import cpu_count
from sys import stdout
from copy import deepcopy
from subprocess import call
from cobra.util import solver
import platform
import pandas as pd
from unipath import Path
from cobra.manipulation.delete import *



# User defined arguments
parser = argparse.ArgumentParser(description='Generate genome-scale metabolic network reconstruction from KEGG BLAST hits.')
parser.add_argument('--input_file', default='none')
parser.add_argument('--file_type', default=7, help='Input file type: fasta=1, diamond blastp output=2, genre sbml=3,clean=5, DeepECv2_result = 6,clean_all = 7')
parser.add_argument('--cleanfile', default='none')
parser.add_argument('--reward', default=0.1, help='reward for new reactions')
parser.add_argument('--iter', default = 3, help='do you want to perform the test suite?')
parser.add_argument('--threshold', default = 8, help='do you want to perform the test suite?')
parser.add_argument('--block_flage', default = 0, help='decrease block reactions: 1, ignore block reactions: 0')
parser.add_argument('--flux_flage', default = 0, help='reward changed by flux of the reactions: 1, ignore reactions fluxes: 0')

parser.add_argument('--media', default='rich', help='List of metabolites composing the media condition. Not required.')
parser.add_argument('--tasks', default=[], help='List of metabolic tasks. Not required.')
parser.add_argument('--org', default='default', help='KEGG organism code. Not required.')
parser.add_argument('--min_frac', default=0.01, help='Minimum objective fraction required during gapfilling')
parser.add_argument('--max_frac', default=0.5, help='Maximum objective fraction allowed during gapfilling')
parser.add_argument('--gram', default='none', help='Type of Gram classificiation (positive or negative)')
parser.add_argument('--out', default='default', help='Name of output GENRE file')
parser.add_argument('--name', default='default', help='ID of output GENRE')
parser.add_argument('--cpu', default=1, help='Number of processors to use')
parser.add_argument('--gapfill', default='yes', help='gapfill your model?')
parser.add_argument('--exchange', default = 1, help='open exchange: 1, shut down exchange: 0')
parser.add_argument('--test', default = 'no', help='do you want to perform the test suite?')

args = parser.parse_args()
with open('/home/kexin/code/bigg/data/biggr2ec.pkl', 'rb') as f:
    biggr2ec = pickle.load(f)
with open('/home/kexin/code/bigg/data/biggec2r.pkl', 'rb') as f:
    biggec2r = pickle.load(f)
print("biggr2ec dictionary loaded successfully.")
with open('/home/kexin/code/bigg/data/seedr2ec.pkl', 'rb') as f:
    seedr2ec = pickle.load(f)
with open('/home/kexin/code/bigg/data/seedec2r.pkl', 'rb') as f:
    seedec2r = pickle.load(f)
print("seedr2ec dictionary loaded successfully.")
# with open('biggr2ec.pkl', 'rb') as f:
#     biggr2ec = pickle.load(f)
# with open('biggec2r.pkl', 'rb') as f:
#     biggec2r = pickle.load(f)
# print("biggr2ec dictionary loaded successfully.")
# with open('seedr2ec.pkl', 'rb') as f:
#     seedr2ec = pickle.load(f)
# with open('seedec2r.pkl', 'rb') as f:
#     seedec2r = pickle.load(f)
# print("seedr2ec dictionary loaded successfully.")

def clean_read_blast(blast_hits):
    pr_d ={}
    with open(blast_hits, 'r') as inFile:
        for line in inFile:
            line = line.split()
            pr_d[line[0]] = line[1]
    return pr_d

def read_clean(input_file):
    threrhold = 0.8
    print('RRRReverse threrhold-->',threrhold)
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
                    # if dis >= threrhold:
                    if dis <= threrhold:
                        try:
                            pr2ec[pr].append(ecid)
                        except KeyError:
                            pr2ec[pr] = [ecid]     
    print('pr2ec-protein number->',len(list(pr2ec.keys())))   
    
    '''print(pr2ec)
    {'fig|837.83.peg.2|PGA7_00000020': ['2.3.1.201'], 'fig|837.83.peg.4|PGA7_00000040': ['2.3.1.286'], 
    'fig|837.83.peg.11|PGA7_00000130': ['4.1.1.81'], 'fig|837.83.peg.12|PGA7_00000140': ['3.1.3.71'], 
    'fig|837.83.peg.14|PGA7_00000160': ['2.7.13.3'], 'fig|837.83.peg.26|PGA7_00000270': ['4.6.1.12'], 
    'fig|837.83.peg.41|PGA7_00000410': ['2.1.2.1'], 'fig|837.83.peg.45|PGA7_00000450': ['2.7.7.41'], 
    'fig|837.83.peg.52|PGA7_00000510': ['2.7.13.3'], 'fig|837.83.peg.57|PGA7_00000550': ['6.3.4.21']}
    '''
    return pr2ec
    
def read_clean_withscore(input_file,threshold=0.5):
    print('threrhold-->',threshold)
    pr2ec = {}
    allec2pr = {}
    ec2pr = {}
    allpr2ec = {}
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
                    ecid = ec.split(':')[-1]
                    # ecid = ec
                    dis = float(dis)
                    try:
                        allpr2ec[pr].update({ecid:dis})
                    except:
                        allpr2ec[pr] = {ecid:dis}
                    try:
                        allec2pr[ecid].update({pr:dis})
                    except:
                        allec2pr[ecid] = {pr:dis}

                    if dis <= threshold:
                        try:
                            ec2pr[ecid].update({pr:dis})
                        except:
                            ec2pr[ecid] = {pr:dis}
                        try:
                            pr2ec[pr].append(ecid)
                        except:
                            pr2ec[pr] = [ecid]
    print('pr2ec-protein number->',len(list(pr2ec.keys())))   
    return pr2ec,ec2pr,allpr2ec,allec2pr

def read_cleandf_withscore(input_file,threshold=8):
    df = pd.read_pickle(input_file)

    allpr2ec = df.to_dict()
    allec2pr = df.T.to_dict()
    
    df = df[df <= threshold]
    # df_cleaned = df.dropna(how='all', axis=0).dropna(how='all', axis=1)

    pr2ec = {col: {row: value for row, value in df[col].dropna().items()} 
                                  for col in df.columns
                                  if not df[col].dropna().empty
                                  }
    ec2pr = {row: {col: value for col, value in df.loc[row].dropna().items()} 
                        for row in df.index 
                        if not df.loc[row].dropna().empty}
    
    return pr2ec,ec2pr,allpr2ec,allec2pr

def maximum_separation(dist_lst, first_grad, use_max_grad):
    opt = 0 if first_grad else -1
    gamma = np.append(dist_lst[1:], np.repeat(dist_lst[-1], 10))
    sep_lst = np.abs(dist_lst - np.mean(gamma))
    sep_grad = np.abs(sep_lst[:-1]-sep_lst[1:])
    if use_max_grad:
        # max separation index determined by largest grad
        max_sep_i = np.argmax(sep_grad)
    else:
        # max separation index determined by first or the last grad
        large_grads = np.where(sep_grad > np.mean(sep_grad))
        max_sep_i = large_grads[-1][opt]
    # if no large grad is found, just call first EC
    if max_sep_i >= 5:
        max_sep_i = 0
    return max_sep_i
'''
for col in df.columns:
    ec = []
    smallest_10_dist_df = df[col].nsmallest(10)
    dist_lst = list(smallest_10_dist_df)
    max_sep_i = maximum_separation(dist_lst, first_grad, use_max_grad)
    for i in range(max_sep_i+1):
        EC_i = smallest_10_dist_df.index[i]
        dist_i = smallest_10_dist_df[i]
'''
def clean2seedr(allpr2ec,threshold):
    universal_scoredict={}
    rxns = {}
    for pr, ec_score in allpr2ec.items():
        smallest_10_dist_df = pd.Series(ec_score).nsmallest(20)
        dist_lst = list(smallest_10_dist_df)
        max_sep_i = maximum_separation(dist_lst, True, True)
        max_sep_ec = [smallest_10_dist_df.index[i] for i in range(max_sep_i+1)]
        
        for ec , score in ec_score.items():
            ec = ec.split(':')[-1]
            if ec in seedec2r.keys():
                rids = seedec2r[ec]
                # print('ec',ec,'|','rids:',rids)
                for rid in rids:
                    rid = rid + '_c'

                    if rid in universal_scoredict.keys():
                        ori = universal_scoredict[rid]
                        universal_scoredict[rid] = min(score,ori)
                        if universal_scoredict[rid] <= threshold:
                            if ori ==  universal_scoredict[rid]:
                                continue
                            elif ec in max_sep_ec:
                                try:
                                    rxns[rid].append(pr) 
                                except KeyError:
                                    rxns[rid] = [pr]
                    else:
                        universal_scoredict[rid] = score
                        if score <= threshold and ec in max_sep_ec:
                            try:
                                rxns[rid].append(pr) 
                            except KeyError:
                                rxns[rid] = [pr]
    print('rxns:',len(rxns))
    # print('rxns:',rxns)
    return universal_scoredict,rxns

def update_predscore(newreactions,allec2pr,newallpr2ec,updateprs,reward,threshold):
    print('len of newreactions:',len(newreactions))
    print('len of updateprs:',len(updateprs))
    # print('before updataer score',[newpredscore[pr] for pr in updateprs])
    count = 0
    for r in newreactions:
        r = r.split('_')[0]
        try:
            ecid = seedr2ec[r]
        except KeyError:
            count+=1
            continue
        for ec in ecid:
            try:
                prd = allec2pr[ec]

            except KeyError:
                count+=1
                continue

            # take the pr with max val
            try:
                pr = min(prd, key=prd.get)
            except ValueError:
                count+=1
                continue

            try:
                s = max(float(newallpr2ec[pr][ec]) - reward, 0)
                if s <= threshold:
                    updateprs.append(pr)

                print(pr,ec,newallpr2ec[pr][ec],'->', s,flush=True)
                newallpr2ec[pr].update({ec: s})
                allec2pr[ec].update({pr: s})

            except KeyError:
                count+=1
    print(f'during update: {count} missing')    
       
    return newallpr2ec,updateprs,allec2pr


def update_predscore_flux(rxn_flux,allec2pr,newpredscore,updateprs,reward,threshold,fluxflage):
    count = 0
    ###
    if fluxflage == 1:
        print('gradient reward by flux...')
    ####
    else:
        print('reward is constant val not effect by flux...',flush=True)

    ####
    print('fluxflage:',fluxflage)
    print('update flux reactions...',flush=True)
    for r in rxn_flux.keys():
        flux = media_flux[r]
        r = r.split('_')[0]
        if fluxflage==1:
            if flux > 0:
                reward = reward*10
            elif flux > 1e-5:
                reward = reward*5
            elif flux > 1e-6:
                reward = reward
        try:
            ecid = seedr2ec[r]
        except KeyError:
            count+=1
            continue
        for ec in ecid:
            try:
                prd = allec2pr[ec]

            except KeyError:
                count+=1
                continue

            # take the pr with max val
            try:
                pr = min(prd, key=prd.get)
            except ValueError:
                count+=1
                continue

            try:
                s = max(float(newpredscore[pr][ec]) - reward, 0)
                if s <= threshold:
                    updateprs.append(pr)
                print(pr,ec,newpredscore[pr][ec],'->', s,flush=True)
                newpredscore[pr].update({ec: s})
                allec2pr[ec].update({pr: s})

            except KeyError:
                count+=1
    print(f'during update: {count} missing')    
       
    return newpredscore,updateprs,allec2pr

def update_predscore_block(blockrxn,allec2pr,newpredscore,updateprs,reward,threshold):
    print('updata block reactions...',flush=True)
    count = 0
    for r in blockrxn:
        r = r.split('_')[0]
        try:
            ecid = seedr2ec[r]
        except KeyError:
            count+=1
            continue
        for ec in ecid:
            try:
                prd = allec2pr[ec]
            except KeyError:
                count+=1
                continue
            # take the pr with max val
            try:
                pr = min(prd, key=prd.get)
            except ValueError:
                count+=1
                continue
            try:
                s = float(newpredscore[pr][ec]) + reward
                # if s <= threshold:
                    # updateprs.append(pr)
                print(pr,ec,newpredscore[pr][ec],'->', s,flush=True)
                newpredscore[pr].update({ec: s})
                allec2pr[ec].update({pr: s})

            except KeyError:
                count+=1
    print(f'during update: {count} missing')    
       
    return newpredscore,updateprs,allec2pr
     
def read_ecpred(input_file):
    pr2ec = {}
    with open(input_file, 'r') as inFile:
        for line in inFile:
            line = line.split()
            pr = line[0]
            items = line[-1].split(';')
            for item in items:
                if item.startswith('None'):
                    break
                elif item.startswith('EC:'):
                    # ec = item.split('/')
                    ecid = item.split(':')[-1]
                    
                    try:
                        pr2ec[pr].append(ecid)
                    except KeyError:
                        pr2ec[pr] = [ecid]     
    print('pr2ec-protein number->',len(list(pr2ec.keys()))) 
    return pr2ec
    
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

# Translate genes to ModelSEED reactions
def _genes_to_rxns(kegg_hits, gene_modelseed, organism):
    ''' Translates genes to ModelSEED reactions '''
    if organism != 'default':
        new_hits = _get_org_rxns(gene_modelseed, organism)
        gene_count = len(kegg_hits)
        kegg_hits |= new_hits
        gene_count = len(kegg_hits) - gene_count
        print('Added', gene_count, 'genes from', organism)

    rxn_db = {}
    for gene in kegg_hits:
        try:
            rxns = gene_modelseed[gene]
        except KeyError:
            continue

        for rxn in rxns:
            rxn = rxn + '_c'
            try:
                rxn_db[rxn].append(gene)
            except KeyError:
                rxn_db[rxn] = [gene]

    return rxn_db

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
    thred = 0.8
    for r in universal_scoredict.keys():
        if len(universal_scoredict[r]) >= thred:
            rxns.append(r)

    return universal_scoredict,rxns

# Get genes for organism from reference genome
def _get_org_rxns(gene_modelseed, organism):
    ''' Get genes for organism from reference genome '''
    rxn_db = {}

    org_genes = []
    for gene in gene_modelseed.keys():
        current = gene.split(':')[0]
        if current == organism:
            org_genes.append(gene)

    return set(org_genes)

# import time.time
# Create draft GENRE and integrate GPRs
def _create_model(rxn_db, universal, input_id):
    ''' Create draft GENRE and integrate GPRs '''
    new_model = cobra.Model('new_model')
    c = 0
    # reactions = []
    tmpuniversal = deepcopy(universal)
    for x in rxn_db.keys():
        # orix = x
        # x = x + '_c'
        # print(orix,'x:',x,'->',rxn_db[orix])
        try:
            rxn = tmpuniversal.reactions.get_by_id(x)
            # rxn = deepcopy(rxn)
            rxn.gene_reaction_rule = ' or '.join(rxn_db[x])
            new_model.add_reactions([rxn])
            c+=1
        except KeyError:
            # print('keyerror:',x)
            continue
        # if c%50 == 0:
        #     # print('time:',time.time(),flush=True)
        #     print('c:',c,flush=True)
    if input_id != 'default':
        new_model.id = input_id
    print('count of reactions:',c)
    return new_model


# Add gene names
def _add_names(model, gene_db):
    ''' Add gene names '''
    for gene in model.genes:
        try:
            gene.name = gene_db[gene.id].title()
        except KeyError:
            continue

    return model


# pFBA gapfiller
def _find_reactions(model, reaction_bag, tasks, obj, fraction, max_fraction, step, file_type):
    ''' pFBA gapfiller that modifies universal reaction bag, removes overlapping reacitons from universal reaction bag
    and resets objective if needed, adds model reaction to universal bag, sets lower bound for metabolic tasks, 
    sets minimum lower bound for previous objective, assemble forward and reverse components of all reactions,
    create objective, based on pFBA, run FBA and identify reactions from universal that are now active'''
    stdout.write('\r[                                         ]')
    stdout.flush()

    # Modify universal reaction bag
    new_rxn_ids = set() #make empty set we will add new reaction ids to
    with reaction_bag as universal: #set the reaction bag as the universal reaction databse

        # Remove overlapping reactions from universal bag, and reset objective if needed
        warnings.filterwarnings('ignore')
        orig_rxn_ids = set()    #original reaction ids start as an empty set
        remove_rxns = []    #reactions to remove is empty vector
        for rxn in model.reactions: #for a reaction in the draft model reactions
            if rxn.id == obj and file_type != 3: #if a reaction is part of the objective function 
                continue

            orig_rxn_ids |= set([rxn.id])
            try:
                test = universal.reactions.get_by_id(rxn.id)
                remove_rxns.append(rxn.id)
            except:
                continue

        # Add model reactions to universal bag
        universal.remove_reactions(list(set(remove_rxns)))
        add_rxns = []
        for x in model.reactions:
            if x.id != obj or file_type == 3:
                add_rxns.append(deepcopy(x))
        universal.add_reactions(add_rxns)

        # Set lower bounds for metaboloic tasks
        if len(tasks) != 0:
            for rxn in tasks:
                try:
                    universal.reactions.get_by_id(rxn).lower_bound = fraction
                except:
                    continue

        stdout.write('\r[---------------                          ]')
        stdout.flush()

        # Set minimum lower bound for previous objective
        universal.objective = universal.reactions.get_by_id(obj) 
        prev_obj_val = universal.slim_optimize()
        print('\nuniversal.objective:', universal.objective)
        print('universal.slim_optimize:', prev_obj_val,flush=True)
        if step == 1:
            print('1 current constraint',universal.reactions.get_by_id(obj))
            prev_obj_constraint = universal.problem.Constraint(universal.reactions.get_by_id(obj).flux_expression, 
        	   lb=prev_obj_val*fraction, ub=prev_obj_val*max_fraction)
        elif step == 2:
            print('2 current constraint',universal.reactions.get_by_id(obj))
            
            prev_obj_constraint = universal.problem.Constraint(universal.reactions.get_by_id(obj).flux_expression, 
               lb=prev_obj_val*max_fraction, ub=prev_obj_val)
        universal.solver.add(prev_obj_constraint)
        universal.solver.update()

        # Assemble forward and reverse components of all reactions
        coefficientDict = {}
        pfba_expr = symengine.RealDouble(0)
        for rxn in universal.reactions:
            if rxn.id in orig_rxn_ids:
                coefficientDict[rxn.forward_variable] = 0.0
                coefficientDict[rxn.reverse_variable] = 0.0
            else:
                coefficientDict[rxn.forward_variable] = 1.0
                coefficientDict[rxn.reverse_variable] = 1.0

        stdout.write('\r[--------------------------               ]')
        stdout.flush()

        # Create objective, based on pFBA
        universal.objective = 0
        universal.solver.update()
        universal.objective = universal.problem.Objective(symengine.RealDouble(0), direction='min', sloppy=True)
        universal.objective.set_linear_coefficients(coefficientDict)
        
        stdout.write('\r[----------------------------------       ]')
        stdout.flush()

        # Run FBA and identify reactions from universal that are now active
        solution = universal.optimize()
        
    # new_rxn_ids = set([rxn.id for rxn in reaction_bag.reactions if abs(solution.fluxes[rxn.id]) > 1e-6]).difference(orig_rxn_ids)
    print('step {}-'.format(step),solution.fluxes[rxn.id])
    print('step {}-'.format(step),len(orig_rxn_ids),list(orig_rxn_ids)[:10])
    coutn=0
    for rxn in reaction_bag.reactions:
        a = solution.fluxes[rxn.id]
        if abs(a) > 1e-13:
            print(rxn.id,'---->',a,end='\n')
        else:
            if a == 0.0:
                coutn+=1
            
    print('step {}-'.format(step),len(reaction_bag.reactions))
    print('step {}-'.format(step),coutn)
    threshold1 = 1e-10
    new_rxn_ids = set([rxn.id for rxn in reaction_bag.reactions if abs(solution.fluxes[rxn.id]) > threshold1]).difference(orig_rxn_ids)
    print('threshold:',threshold1)
    stdout.write('\r[-----------------------------------------]\n')
    warnings.filterwarnings('default')
    print('step {} --new_rxn_ids:==>'.format(step),len(new_rxn_ids),'\n',new_rxn_ids)
    return(new_rxn_ids)     


# pFBA gapfiller
def weighted_find_reactions(model, reaction_bag, tasks, obj, fraction, max_fraction, step, file_type,weight_dict):
    ''' pFBA gapfiller that modifies universal reaction bag, removes overlapping reacitons from universal reaction bag
    and resets objective if needed, adds model reaction to universal bag, sets lower bound for metabolic tasks, 
    sets minimum lower bound for previous objective, assemble forward and reverse components of all reactions,
    create objective, based on pFBA, run FBA and identify reactions from universal that are now active'''
    # stdout.write('\r[                                         ]')
    # stdout.flush()

    # Modify universal reaction bag
    new_rxn_ids = set() #make empty set we will add new reaction ids to
    with reaction_bag as universal: #set the reaction bag as the universal reaction databse

        # Remove overlapping reactions from universal bag, and reset objective if needed
        print('#Remove overlapping reactions from universal bag, and reset objective if needed')
        warnings.filterwarnings('ignore')
        orig_rxn_ids = set()    #original reaction ids start as an empty set
        remove_rxns = []    #reactions to remove is empty vector
        for rxn in model.reactions: #for a reaction in the draft model reactions
            if rxn.id == obj and file_type != 3: #if a reaction is part of the objective function 
                continue

            orig_rxn_ids |= set([rxn.id])
            try:
                test = universal.reactions.get_by_id(rxn.id)
                remove_rxns.append(rxn.id)
            except:
                continue

        # Add model reactions to universal bag
        print('#add model reactions to universal bag')
        universal.remove_reactions(list(set(remove_rxns)))
        add_rxns = []
        for x in model.reactions:
            if x.id != obj or file_type == 3:

                add_rxns.append(x.copy())
        universal.add_reactions(add_rxns)
      


        # Set lower bounds for metaboloic tasks
        print('#set lower bounds for metaboloic tasks')
        if len(tasks) != 0:
            for rxn in tasks:
                try:
                    universal.reactions.get_by_id(rxn).lower_bound = fraction
                except:
                    continue

        # stdout.write('\r[---------------                          ]')
        # stdout.flush()

        # Set minimum lower bound for previous objective
        print('#set minimum lower bound for previous objective')
        universal.objective = universal.reactions.get_by_id(obj) 
        prev_obj_val = universal.slim_optimize()
        print('#prev_obj_val:',prev_obj_val)
        # if prev_obj_val < 0:
        #     prev_obj_val = abs(prev_obj_val)
        if step == 1:
            prev_obj_constraint = universal.problem.Constraint(universal.reactions.get_by_id(obj).flux_expression, 
        	   lb=prev_obj_val*fraction, ub=prev_obj_val*max_fraction)
        elif step == 2:
            prev_obj_constraint = universal.problem.Constraint(universal.reactions.get_by_id(obj).flux_expression, 
               lb=prev_obj_val*max_fraction, ub=prev_obj_val)
        universal.solver.add(prev_obj_constraint)
        universal.solver.update()

        # Assemble forward and reverse components of all reactions
        print('#assemble forward and reverse components of all reactions')
        coefficientDict = {}


        pfba_expr = symengine.RealDouble(0)
        for rxn in universal.reactions:
            if rxn.id in orig_rxn_ids:
                coefficientDict[rxn.forward_variable] = 0.005
                coefficientDict[rxn.reverse_variable] = 0.005
            else:
                if rxn.id in weight_dict:
                    # print('rxn.id:',rxn.id,'weight_dict[rxn.id]:',weight_dict[rxn.id],rxn.id in rm_rxn)
                    # coefficientDict[rxn.forward_variable] = weight_dict[rxn.id]
                    # coefficientDict[rxn.reverse_variable] = weight_dict[rxn.id]
                    # coefficientDict[rxn.forward_variable] = 90 * (1 - (weight_dict[rxn.id] - min_old) / (max_old - min_old))
                    # coefficientDict[rxn.reverse_variable] = 90 * (1 - (weight_dict[rxn.id] - min_old) / (max_old - min_old))
                    # coefficientDict[rxn.forward_variable] = 1 * 1/weight_dict[rxn.id]
                    # coefficientDict[rxn.reverse_variable] = 1 * 1/weight_dict[rxn.id]
                    # print('rxn.id:',rxn.id,'weight_dict[rxn.id]:',weight_dict[rxn.id])
                    # coefficientDict[rxn.forward_variable] = 1 - weight_dict[rxn.id]
                    # coefficientDict[rxn.reverse_variable] = 1 - weight_dict[rxn.id]
                    # RRREVERSED
                    coefficientDict[rxn.forward_variable] = weight_dict[rxn.id]
                    coefficientDict[rxn.reverse_variable] = weight_dict[rxn.id]
                    ##normalize penalty score
                    # coefficientDict[rxn.forward_variable] = weight_dict[rxn.id]*100
                    # coefficientDict[rxn.reverse_variable] = weight_dict[rxn.id]
                    # gapseq = ((100-0)/(1-0.3)) * (weight_dict[rxn.id] - 0.3)
                else:
                    coefficientDict[rxn.forward_variable] = 100
                    coefficientDict[rxn.reverse_variable] = 100
                # else:##normalize penalty score
                #     coefficientDict[rxn.forward_variable] = 100
                #     coefficientDict[rxn.reverse_variable] = 100
        # sss =  
        # print('score:',sss)
        # effrs = set()
        # for orir in orig_rxn_ids:
        #     ec = biggr2ec[orir]
        #     for curr_ec in ec:
        #         havers = biggec2r[curr_ec]
        #         for nrs in havers:
        #             # remove nrs in weight_dict
        #             if nrs in weight_dict:
        #                 effrs.add(nrs)
        #                 # try:
        #                 #     coefficientDict[universal.reactions.get_by_id(nrs).forward_variable] += 2
        #                 #     coefficientDict[universal.reactions.get_by_id(nrs).reverse_variable] += 2
        #                 # except:
        #                 #     continue
        # # print('effrs:',len(effrs),effrs)
        # for effr in effrs:
        #     if effr in rm_rxn:
        #         print('effr:',effr)
        #     try:
        #         coefficientDict[universal.reactions.get_by_id(effr).forward_variable] += 10
        #         coefficientDict[universal.reactions.get_by_id(effr).reverse_variable] += 10
        #     except:
        #         continue 
  

        # Create objective, based on pFBA
        print('#create objective, based on pFBA')
        universal.objective = 0
        universal.solver.update()
        universal.objective = universal.problem.Objective(symengine.RealDouble(0), direction='min', sloppy=True)
        universal.objective.set_linear_coefficients(coefficientDict)
        
        # stdout.write('\r[----------------------------------       ]')
        # stdout.flush()

        # Run FBA and identify reactions from universal that are now active
        print('#run FBA and identify reactions from universal that are now active')
        solution = universal.optimize()
        print('# run FBA solution:',solution)
    # coutn=0
    # for rxn in reaction_bag.reactions:
    #     a = solution.fluxes[rxn.id]
    #     # if rxn.id in orig_rxn_ids:
    #     #     continue    
    #     if abs(a) > 1e-13 and rxn.id not in orig_rxn_ids:
    #         print(rxn.id,'---->',a,end='\n')
    #     else:
    #         if a == 0.0:
    #             coutn+=1
    # print('count',coutn)
    fluxthreshold = 1e-6
    print('fluxthreshold:',fluxthreshold)
    # new_rxn_ids = set([rxn.id for rxn in reaction_bag.reactions if abs(solution.fluxes[rxn.id]) > 1e-6]).difference(orig_rxn_ids)
    # new_rxn_ids = set([rxn.id for rxn in reaction_bag.reactions if abs(solution.fluxes[rxn.id]) > fluxthreshold]).difference(orig_rxn_ids)
    new_rxn_flux = {}
    new_rxn_ids = set()
    for rxn in reaction_bag.reactions:
        if abs(solution.fluxes[rxn.id]) > fluxthreshold:
            if rxn.id not in orig_rxn_ids:
                new_rxn_ids.add(rxn.id)
                new_rxn_flux[rxn.id] = solution.fluxes[rxn.id]
    # stdout.write('\r[-----------------------------------------]\n')
    print('len(new_rxn_ids)',len(new_rxn_ids),'\n',new_rxn_ids,flush=True)
    warnings.filterwarnings('default')
    print('return with flux:',len(new_rxn_flux.keys()))
    return new_rxn_ids ,coefficientDict,new_rxn_flux


# Add new reactions to model
def _gapfill_model(model, universal, new_rxn_ids, obj, step):
    '''Adds new reactions to model by getting reactions and metabolites to be added to the model, creates gapfilled model, 
    and identifies extracellular metabolites that still need exchanges '''
    # Get reactions and metabolites to be added to the model
    new_rxns = []
    if step == 1: new_rxns.append(deepcopy(universal.reactions.get_by_id(obj)))
    for rxn in new_rxn_ids: 
        if rxn != obj:
            new_rxns.append(deepcopy(universal.reactions.get_by_id(rxn)))
    
    # Create gapfilled model 
    model.add_reactions(new_rxns)
    model.objective = model.problem.Objective(model.reactions.get_by_id(obj).flux_expression, direction='max')

    # Identify extracellular metabolites still need exchanges
    for cpd in model.metabolites:
        if cpd.compartment != 'extracellular':
            continue
        else:
            try:
                test = model.reactions.get_by_id('EX_' + cpd.id)
            except KeyError:
                exch_id = 'EX_' + cpd.id
                model.add_boundary(cpd, type='exchange', reaction_id=exch_id, lb=-1000.0, ub=1000.0)
                model.reactions.get_by_id(exch_id).name = cpd.name + ' exchange'

    return model


# Set uptake of specific metabolites in complete medium gap-filling
def _set_base_inputs(model, universal):
    tasks = ['EX_cpd00035_e','EX_cpd00051_e','EX_cpd00132_e','EX_cpd00041_e','EX_cpd00084_e','EX_cpd00053_e','EX_cpd00023_e',
    'EX_cpd00033_e','EX_cpd00119_e','EX_cpd00322_e','EX_cpd00107_e','EX_cpd00039_e','EX_cpd00060_e','EX_cpd00066_e','EX_cpd00129_e',
    'EX_cpd00054_e','EX_cpd00161_e','EX_cpd00065_e','EX_cpd00069_e','EX_cpd00156_e','EX_cpd00027_e','EX_cpd00149_e','EX_cpd00030_e',
    'EX_cpd00254_e','EX_cpd00971_e','EX_cpd00063_e','EX_cpd10515_e','EX_cpd00205_e','EX_cpd00099_e']

    new_rxns = []
    for exch in tasks: 
        try:
            test = model.reactions.get_by_id(exch)
        except:
            new_rxns.append(deepcopy(universal.reactions.get_by_id(exch)))
    model.add_reactions(new_rxns)
    for exch in tasks: model.reactions.get_by_id(exch).bounds = (-1000., -0.01)

    return model


def _add_annotation(model, gram, obj='built'):
    ''' Add gene, metabolite, reaction ,biomass reaction annotations '''
    # Genes
    for gene in model.genes:
        gene._annotation = {}
        gene.annotation['sbo'] = 'SBO:0000243'
        gene.annotation['kegg.genes'] = gene.id
    
    # Metabolites
    for cpd in model.metabolites: 
        cpd._annotation = {}
        cpd.annotation['sbo'] = 'SBO:0000247'
        if 'cpd' in cpd.id: cpd.annotation['seed.compound'] = cpd.id.split('_')[0]

    # Reactions
    for rxn in model.reactions:
        rxn._annotation = {}
        if 'rxn' in rxn.id: rxn.annotation['seed.reaction'] = rxn.id.split('_')[0]
        compartments = set([x.compartment for x in list(rxn.metabolites)])
        if len(list(rxn.metabolites)) == 1:
            rxn.annotation['sbo'] = 'SBO:0000627' # exchange
        elif len(compartments) > 1:
            rxn.annotation['sbo'] = 'SBO:0000185' # transport
        else:
            rxn.annotation['sbo'] = 'SBO:0000176' # metabolic

    # Biomass reactions
    if obj == 'built':
        try:
            model.reactions.EX_biomass.annotation['sbo'] = 'SBO:0000632'
        except:
            pass
        if gram == 'none':
            biomass_ids = ['dna_rxn','rna_rxn','protein_rxn','teichoicacid_rxn','lipid_rxn','cofactor_rxn','rxn10088_c','biomass_rxn']
        else:
            biomass_ids = ['dna_rxn','rna_rxn','protein_rxn','teichoicacid_rxn','peptidoglycan_rxn','lipid_rxn','cofactor_rxn','GmPos_cellwall','rxn10088_c','GmNeg_cellwall','biomass_rxn_gp','biomass_rxn_gn']
        for x in biomass_ids:
            try:
                model.reactions.get_by_id(x).annotation['sbo'] = 'SBO:0000629'
            except:
                continue
    else:
        model.reactions.get_by_id(obj).annotation['sbo'] = 'SBO:0000629'

    return model
    

# Run some basic checks on new models
def _checkModel(pre_reactions, pre_metabolites, post_model):
    print('\n\tChecking new model...',flush=True)  
    ''' Run basic checks on new models (checking for objective flux'''

	# Check for objective flux
    new_genes = len(post_model.genes)
    new_rxn_ids = set([x.id for x in post_model.reactions]).difference(pre_reactions)
    new_cpd_ids = set([x.id for x in post_model.metabolites]).difference(pre_metabolites)
    test_flux = round(post_model.slim_optimize(), 3)

	
	# Report to user
    print('\tDraft reconstruction had', str(new_genes), 'genes,', str(len(pre_reactions)), 'reactions, and', str(len(pre_metabolites)), 'metabolites')
    print('\tGapfilled', str(len(new_rxn_ids)), 'reactions and', str(len(new_cpd_ids)), 'metabolites\n')
    print('\tFinal reconstruction has', str(len(post_model.reactions)), 'reactions and', str(len(post_model.metabolites)), 'metabolites')
    print('\tFinal objective flux is', str(round(test_flux, 3)))

    return test_flux


#----------------------------------------------------------------------------------------------------------------------#
if __name__ == "__main__":
    # Process input settings
    input_file = str(args.input_file)
    cleanfile = str(args.cleanfile)
    out_file = str(args.out)
    file_type = int(args.file_type)
    name = str(args.name)
    org = str(args.org)
    inter = int(args.iter)
    reward = float(args.reward)
    threshold = float(args.threshold)
    block_flage = int(args.block_flage)
    flux_flage = int(args.flux_flage)
    try:
        media = str(args.media).split(",")
    except:
        media = str(args.media)
    min_frac = float(args.min_frac)
    max_frac = float(args.max_frac)
    metabolic_tasks = list(args.tasks)
    new_id = str(args.name)
    gram_type = str(args.gram)
    processors = int(args.cpu)
    gapfill = str(args.gapfill)
    exchange_arg = int(args.exchange)
    test = str(args.test)
    print('all args:',args,flush=True)
    if gram_type == 'positive':
        print('\nUsing Gram positive objective function')
        universal_obj = 'biomass_GmPos'
    elif gram_type == 'negative':
        print('\nUsing Gram negative objective function')
        universal_obj = 'biomass_GmNeg'
    else:
        universal_obj = 'biomass'

    if min_frac <= 0.0 or min_frac > 1.0:
        print('WARNING: Improper minimum fraction selected. Defaulting to 0.01')
        min_frac = 0.01
    if max_frac <= 0.0 or max_frac > 1.0:
        print('WARNING: Improper maximum fraction selected. Defaulting to 0.5')
        max_frac = 0.5
    if max_frac < min_frac:
        print('WARNING: Input maximum fraction less than minimum fraction. Minimum set to half maximum')
        min_frac = max_frac * 0.5

    if org != 'default':
        print('Including additional genes from KEGG genome of', org)

    # Maximum fraction should not be too high, otherwise the gapfiller adds too many reactions
    print('Using minimum objective flux fraction of', min_frac,'and maximum fraction of', max_frac)

    if processors > cpu_count():
        print('WARNING: Requested more processors than are available. Using maximum of', cpu_count())
        processors = cpu_count()
    print('Using', processors, 'processor(s)\n')


    # Load databases
    script_path = '/ibex/user/niuk0a/anaconda3/envs/recon/lib/python3.9/site-packages/reconstructor'
    # /home/kexin/anaconda3/envs/recon
    script_path = '/home/kexin/anaconda3/envs/recon/lib/python3.9/site-packages/reconstructor'
    print('Loading GENRE construction databases...')
    # script_path = str(os.path.dirname(os.path.realpath(__file__)))
    kegg_prot_db = script_path + '/refs/screened_kegg_prokaryotes_pep_db'
    stdout.write('\r[                                         ]')
    stdout.flush()
    filename = script_path + '/refs/gene_modelseed.pickle'
    with open(filename, 'rb') as f: gene_modelseed = pickle.load(f)
    stdout.write('\r[---------------                          ]')
    stdout.flush()
    filename = script_path + '/refs/universal.pickle'
    with open(filename, 'rb') as f: universal = pickle.load(f)
    stdout.write('\r[------------------------------           ]')
    stdout.flush()
    filename = script_path + '/refs/gene_names.pickle'
    with open(filename, 'rb') as f: gene_names = pickle.load(f)
    # stdout.write('\r[-----------------------------------------]\n')
    print('-->GENRE construction databases loaded\n',flush=True)

    # Check input file type
    # if file_type == 1:
    #     print('Aligning peptide sequences to KEGG database, may take some time...')
    #     # blast_results = input_file.rstrip('fastn') + 'KEGGprot.out'
    #     # print('Saving BLASTp results to', blast_results,'\n')
    #     # _run_blast(input_file, blast_results, kegg_prot_db, str(processors), script_path)
    # elif file_type == 2:
    #     blast_results = input_file
    # elif file_type == 5:
    ## inter

    for i in range(1,int(inter)+1):
        print('*'*50)
        print('Inter:',i)
        print('inter=',inter)
        print('*'*50)
        
        if file_type == 5:
            if i == 1:
                print('Reading clean results from', input_file,flush=True)
                fasta_file = input_file
                clean_file = cleanfile
                # fasta_file = input_file.rstrip('_maxsep.csv') + '.faa'
                pr2gene = {}
                # gene_hits = read_blast_the(blast_results,the=80)
                print('gene_hits-->the = 80')
                # pr2ec = read_clean(input_file)

                start_time = time.time()
                pr2ec,ec2pr,predscore,allec2pr = read_clean_withscore(clean_file,threshold=threshold)
                print('pr2ec number->',len(list(pr2ec.keys())),flush=True)
                print('ec2pr number->',len(list(ec2pr.keys())),flush=True)
                print('predscore number->',len(list(predscore.keys())),flush=True)
                print('time cost for function read_clean_withscore:',time.time()-start_time)
                start_time = time.time()
                universal_scoredict,rxns = clean2seedr(predscore,threshold=threshold)
                print('time cost for function clean2seedr:',time.time()-start_time)
                newpredscore = predscore
                updateprs = []
            else:   
                if len(updateprs) == 0:
                    print('No updateprs')
                    print(f'End No more iteration no more {i}############',flush=True)
                    break
                predscore = newpredscore
                updateprs = []
                universal_scoredict,rxns = clean2seedr(predscore,threshold=threshold)
        if file_type == 7:
            if i == 1:
                print('Reading clean results from', input_file,flush=True)     
                print('NOTICE: current using distance map from clean',flush=True)
                fasta_file = input_file
                clean_file = cleanfile
                
                start_time = time.time()
                pr2ec,ec2pr,predscore,allec2pr = read_cleandf_withscore(clean_file,threshold=threshold)
                print('time cost for function read_clean_withscore:',time.time()-start_time)
                start_time = time.time()
                universal_scoredict,rxns = clean2seedr(predscore,threshold=threshold)
                print('time cost for function clean2seedr:',time.time()-start_time)
                newpredscore = predscore
                updateprs = []
            else:
                if len(updateprs) == 0:
                    print('No updateprs')
                    print(f'End No more iteration no more {i}############',flush=True)
                    break
                predscore = newpredscore
                updateprs = []
                universal_scoredict,rxns = clean2seedr(predscore,threshold=threshold)


        elif file_type == 6:
            fasta_file = input_file.rstrip('.DeepECv2_result.txt') + '.PATRIC.faa'
            pr2gene = {}
            pr2ec = read_ecpred(input_file)
        else:
            try:
                draft_genre = cobra.io.read_sbml_model(input_file)
            except:
                draft_genre = cobra.io.load_json_model(input_file)

        if file_type == 5 or file_type == 6 or file_type == 7:

            print('Draft clean_to_rxns has', len(rxns), 'reactions',flush=True)
            draft_genre = _create_model(rxns, universal, new_id)
            if block_flage == 1:
                block_rxns = cobra.flux_analysis.find_blocked_reactions(draft_genre)
                print('Blocked reactions:', len(block_rxns),flush=True)
                dead_end_rxns = block_rxns
                print("Reactions involving dead-end metabolites:", len(dead_end_rxns))

                newpredscore, updateprs, allec2pr = update_predscore_block(dead_end_rxns,allec2pr,newpredscore,updateprs,reward,threshold)
            print('Draft draft_genre has', len(draft_genre.reactions), 'reactions',flush=True)
            draft_genre = _add_names(draft_genre, gene_names)
            print('Draft draft_genre_add_names has', len(draft_genre.reactions), 'reactions',flush=True)
        else:
            universal_obj = str(draft_genre.objective.expression).split()[0].split('*')[-1]
        

        # Handle media conditions
        if media == 'rich':
            media = ['cpd00001_e','cpd00035_e','cpd00041_e','cpd00023_e','cpd00119_e','cpd00107_e','cpd00060_e','cpd00161_e','cpd00069_e','cpd00084_e','cpd00033_e'
        'cpd00322_e','cpd00066_e','cpd00054_e','cpd00065_e','cpd00156_e','cpd00220_e','cpd00644_e','cpd00393_e','cpd00133_e','cpd00263_e','cpd00104_e','cpd00149_e',
        'cpd00971_e','cpd00099_e','cpd00205_e','cpd00009_e','cpd00063_e','cpd00254_e','cpd10515_e','cpd00030_e','cpd00242_e','cpd00226_e','cpd01242_e','cpd00307_e',
        'cpd00092_e','cpd00117_e','cpd00067_e''cpd00567_e','cpd00132_e','cpd00210_e','cpd00320_e','cpd03279_e','cpd00246_e','cpd00311_e','cpd00367_e','cpd00277_e',
        'cpd00182_e','cpd00654_e','cpd00412_e','cpd00438_e','cpd00274_e','cpd00186_e','cpd00637_e','cpd00105_e','cpd00305_e','cpd00309_e','cpd00098_e','cpd00207_e',
        'cpd00082_e','cpd00129_e']
        elif media == 'minimal':
            media = ['cpd00001_e','cpd00065_e','cpd00060_e','cpd00322_e','cpd00129_e','cpd00156_e','cpd00107_e','cpd00084_e', 
        'cpd00149_e','cpd00099_e','cpd10515_e','cpd00030_e','cpd00254_e','cpd00063_e','cpd00205_e','cpd00009_e','cpd00971_e','cpd00242_e',
        'cpd00104_e','cpd00644_e','cpd00263_e','cpd00082_e']
        elif media == 'default':
            media = ['cpd00035_e','cpd00051_e','cpd00132_e','cpd00041_e','cpd00084_e','cpd00053_e','cpd00023_e',
        'cpd00033_e','cpd00119_e','cpd00322_e','cpd00107_e','cpd00039_e','cpd00060_e','cpd00066_e','cpd00129_e',
        'cpd00054_e','cpd00161_e','cpd00065_e','cpd00069_e','cpd00156_e','cpd00027_e','cpd00149_e','cpd00030_e',
        'cpd00254_e','cpd00971_e','cpd00063_e','cpd10515_e','cpd00205_e','cpd00099_e']
        elif media =='LB':
            media = [
        'cpd00001_e','cpd00007_e','cpd00009_e','cpd00018_e','cpd00023_e','cpd00027_e','cpd00028_e','cpd00030_e','cpd00033_e','cpd00034_e','cpd00035_e',
        'cpd00039_e','cpd00041_e','cpd00046_e','cpd00048_e','cpd00051_e','cpd00054_e','cpd00058_e','cpd00060_e','cpd00063_e','cpd00065_e','cpd00066_e',
        'cpd00067_e','cpd00069_e','cpd00084_e','cpd00091_e','cpd00092_e','cpd00099_e','cpd00107_e','cpd00119_e','cpd00126_e','cpd00129_e','cpd00149_e',
        'cpd00156_e','cpd00161_e','cpd00182_e','cpd00184_e','cpd00205_e','cpd00215_e','cpd00218_e','cpd00219_e','cpd00220_e','cpd00226_e','cpd00239_e',
        'cpd00246_e','cpd00249_e','cpd00254_e','cpd00311_e','cpd00322_e','cpd00381_e','cpd00383_e','cpd00393_e','cpd00438_e','cpd00531_e','cpd00541_e',
        'cpd00644_e','cpd00654_e','cpd00793_e','cpd00971_e','cpd01012_e','cpd01048_e','cpd03424_e','cpd10515_e','cpd10516_e','cpd11595_e'    
            ]
        elif media == 'NB': ## add 'cpd00020_e'(pyruvate) to the media LB 
            media = [
        'cpd00001_e','cpd00007_e','cpd00009_e','cpd00018_e','cpd00023_e','cpd00027_e','cpd00028_e','cpd00030_e','cpd00033_e','cpd00034_e','cpd00035_e',
        'cpd00039_e','cpd00041_e','cpd00046_e','cpd00048_e','cpd00051_e','cpd00054_e','cpd00058_e','cpd00060_e','cpd00063_e','cpd00065_e','cpd00066_e',
        'cpd00067_e','cpd00069_e','cpd00084_e','cpd00091_e','cpd00092_e','cpd00099_e','cpd00107_e','cpd00119_e','cpd00126_e','cpd00129_e','cpd00149_e',
        'cpd00156_e','cpd00161_e','cpd00182_e','cpd00184_e','cpd00205_e','cpd00215_e','cpd00218_e','cpd00219_e','cpd00220_e','cpd00226_e','cpd00239_e',
        'cpd00246_e','cpd00249_e','cpd00254_e','cpd00311_e','cpd00322_e','cpd00381_e','cpd00383_e','cpd00393_e','cpd00438_e','cpd00531_e','cpd00541_e',
        'cpd00644_e','cpd00654_e','cpd00793_e','cpd00971_e','cpd01012_e','cpd01048_e','cpd03424_e','cpd10515_e','cpd10516_e','cpd11595_e','cpd00020_e'   
            ]
        elif media == 'che':
            simulated_media = [
            'cpd00020_e',  # Pyruvate
            'cpd30698_e',  # Casamino acids
            'cpd26232_e',  # calcium chloride
            'cpd09396_e',  # Potassium sulfate
            'cpd00254_e',  # Magnesium
            'cpd09400_e',  # Sodium sulfate
            'cpd20826_e',  # Silica
            'cpd09695_e',  # Sr+
            'cpd09225_e',  # Borate
            'cpd00966_e',  # Bromide
            'cpd19118_e',  # Ammonium chloride
            'cpd00075_e',  # NO2; NO2-; Nitrite; nitrite
            'cpd17321_e',  # Rubidium Rb+; Rubidium; Rubidium cation; Rubidium ion; Rubidium(1+); rubidium ion
            'cpd27384_e',  # Li+; lithium ion
            'cpd00001_e', 'cpd00007_e', 'cpd00009_e', 'cpd00018_e', 'cpd00023_e', 'cpd00027_e', 'cpd00028_e', 'cpd00030_e', 'cpd00033_e', 'cpd00034_e', 'cpd00035_e',
            'cpd00039_e', 'cpd00041_e', 'cpd00046_e', 'cpd00048_e', 'cpd00051_e', 'cpd00054_e', 'cpd00058_e', 'cpd00060_e', 'cpd00063_e', 'cpd00065_e', 'cpd00066_e',
            'cpd00067_e', 'cpd00069_e', 'cpd00084_e', 'cpd00091_e', 'cpd00092_e', 'cpd00099_e', 'cpd00107_e', 'cpd00119_e', 'cpd00126_e', 'cpd00129_e', 'cpd00149_e',
            'cpd00156_e', 'cpd00161_e', 'cpd00182_e', 'cpd00184_e', 'cpd00205_e', 'cpd00215_e', 'cpd00218_e', 'cpd00219_e', 'cpd00220_e', 'cpd00226_e', 'cpd00239_e',
            'cpd00246_e', 'cpd00249_e', 'cpd00254_e', 'cpd00311_e', 'cpd00322_e', 'cpd00381_e', 'cpd00383_e', 'cpd00393_e', 'cpd00438_e', 'cpd00531_e', 'cpd00541_e',
            'cpd00644_e', 'cpd00654_e', 'cpd00793_e', 'cpd00971_e', 'cpd01012_e', 'cpd01048_e', 'cpd03424_e', 'cpd10515_e', 'cpd10516_e', 'cpd11595_e', 'cpd00020_e'
            ]
        else:
            media = media
        print(media)
        
        # Set media condition
        if len(media) != 0:
            media_condition = set(['EX_' + cpd for cpd in media])
            universal_reactions = set([x.id for x in universal.reactions])
            for rxn in universal_reactions:
                if rxn.startswith('EX_') == True:
                    universal.reactions.get_by_id(rxn).bounds = (0, 1000.0)
                if rxn in media_condition:
                    universal.reactions.get_by_id(rxn).bounds = (-1000.0, 10000)

        # Gapfill new model
        if gapfill == 'yes':
            if file_type != 3:
                print('Identifying new metabolism (Step 1 of 2)...')
            if file_type == 3:
                print('Identifying new metabolism...')
                
            draft_reactions = set([x.id for x in draft_genre.reactions])
            draft_metabolites = set([x.id for x in draft_genre.metabolites])
            warnings.filterwarnings('ignore')
            # new_reactions = _find_reactions(draft_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 1, file_type)
            new_reactions ,coefficientDict,new_rxn_flux = weighted_find_reactions(draft_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 1, file_type,universal_scoredict)
            # new_reactions,coefficientDict = weighted_find_reactions(draft_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 1, file_type,universal_scoredict)
            print('new_reactions id',new_reactions)
            # newscoredict,updateprs = update_cleanscore(new_reactions,newscoredict,updateprs)
            #newrflux ?++ for rewards#
        
            # newpredscore, updateprs, allec2pr = update_predscore(new_reactions,allec2pr,newpredscore,updateprs,reward,threshold)
            newpredscore, updateprs, allec2pr = update_predscore_flux(new_rxn_flux,allec2pr,newpredscore,updateprs,reward,threshold,flux_flage)
            print('af update new_reactions',new_reactions)
            print('updateprs:',updateprs)
            ## remove unsatisfactory reactions
            if gram_type == 'positive':
                new_reactions = set([x for x in new_reactions if x not in ['GmNeg_cellwall','biomass_GmNeg']])
            elif gram_type == 'negative':
                new_reactions = set([x for x in new_reactions if x not in ['GmPos_cellwall','biomass_GmPos']])
            print('new_reactions aft gram',new_reactions)
            filled_genre = _gapfill_model(draft_genre, universal, new_reactions, universal_obj, 1)
            
            if file_type != 3:
                print('Identifying new metabolism (Step 2 of 2)...')
                filled_genre = _set_base_inputs(filled_genre, universal)
                # media_reactions,coefficientDict = weighted_find_reactions(filled_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 2, file_type,universal_scoredict)
                media_reactions,coefficientDict,media_flux = weighted_find_reactions(filled_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 2, file_type,universal_scoredict)

                print('bf update media new_reactions',media_reactions)
                
                # media_reactions = _find_reactions(filled_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 2, file_type)
                ## remove unsatisfactory reactions
                ## maps media_reaction with one protein in the clean file
                # newpredscore,updateprs,allec2pr = update_predscore(media_reactions,allec2pr,newpredscore,updateprs,reward,threshold)
                newpredscore,updateprs,allec2pr = update_predscore_flux(media_flux,allec2pr,newpredscore,updateprs,reward,threshold,flux_flage)
                print('aft update media new_reactions',len(media_reactions),media_reactions)
                print('updateprs:',len(updateprs),updateprs)
                
                if gram_type == 'positive':
                    media_reactions = set([x for x in media_reactions if x not in ['GmNeg_cellwall','biomass_GmNeg']])
                elif gram_type == 'negative':
                    media_reactions = set([x for x in media_reactions if x not in ['GmPos_cellwall','biomass_GmPos']])
                final_genre = _gapfill_model(filled_genre, universal, media_reactions, universal_obj, 2)
                final_genre = _add_annotation(final_genre, gram_type)
            else: 
                final_genre = _add_annotation(filled_genre, universal_obj)
        else:
            draft_reactions = set([x.id for x in draft_genre.reactions])
            draft_metabolites = set([x.id for x in draft_genre.metabolites])
            final_genre = draft_genre
            final_genre = _add_annotation(final_genre, gram_type)
        warnings.filterwarnings('default')

        print('# Correct exchanges and check new model')
        if exchange_arg == 0:
            for exch in final_genre.exchanges: exch.bounds = (0., 0.)
        else:
            for exch in final_genre.exchanges: exch.bounds = (-1000., 1000.)
        for rxn in final_genre.reactions:
            if 'Exchange reaction for' in rxn.name:
                rxn.name = list(rxn.metabolites)[0].name + ' exchange'
        biomass = _checkModel(draft_reactions, draft_metabolites, final_genre)
        
        print('# Write new model to sbml')
        input_file = input_file.split('/')[-1] # write to working directory
        if file_type == 1:
            if new_id != 'default':
                out_file = input_file.rstrip('fastn') + new_id + '.sbml'
            else:
                out_file = input_file.rstrip('fastn') + 'sbml'
        elif file_type == 2:
            if new_id != 'default':
                if input_file != 'none':
                    out_file = input_file.rstrip('out') + new_id + '.sbml'
                else:
                    out_file = new_id + '.sbml'
            else:
                if org != 'default':
                    out_file = org + '.sbml'
                else:
                    out_file = input_file.rstrip('out') + 'sbml'
        elif file_type == 3:
            if new_id != 'default':
                out_file = input_file.rstrip('sbml') + new_id + '.extended.sbml'
            else:
                out_file = input_file.rstrip('sbml') + 'extended.sbml'
        
        elif file_type == 5:
            if new_id != 'default':
                out_file = cleanfile.rstrip('.csv') + new_id + f'iter_{i}'+ '.sbml'
            else:
                out_file = cleanfile.rstrip('.csv') + 'test' +  f'iter_{i}'+'.sbml'
                
        elif file_type == 6:
            if new_id != 'default':
                out_file = input_file.rstrip('_result.txt') + new_id + '.sbml'
            else:
                out_file = input_file.rstrip('_result.txt') + 'test' + '.sbml'
        elif file_type == 7:
            if new_id != 'default':
                out_file = cleanfile.rstrip('.pkl') + new_id + f'iter_{i}'+ '.sbml'
            else:
                out_file = cleanfile.rstrip('.pkl') + 'test' +  f'iter_{i}'+'.sbml'

        print('\nSaving new GENRE to', out_file, '\n')
        ## save the dictionary
        thr = str(threshold).split('.')[-1]
        if float(threshold) >= 1.0:
            thr = str(float(threshold)).split('.')[0]
        else:
            thr = str(threshold).split('.')[-1]
        with open(f'{name}_t{thr}_newpredscore_{i}.pkl', 'wb') as f:
            pickle.dump(newpredscore, f)
        print(f'saved newpredscore {i}')
        with open(f'{name}_t{thr}_updateprs_{i}.pkl', 'wb') as f:
            pickle.dump(updateprs, f)
        print(f'saved updateprs {i}')
        cobra.io.write_sbml_model(final_genre, out_file)
        print('*'*50)
        print('End of Inter:',i)
        print('*'*50)
