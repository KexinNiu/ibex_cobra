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
from funcarve_utils import *

# Ensure argparse is imported
import argparse

# Parse arguments
args = argparse.ArgumentParser(description='Generate genome-scale metabolic network reconstruction from KEGG BLAST hits.')
args.add_argument('--input_file', default='none')
args.add_argument('--file_type', default=1, help='Input file type: cleandf=1,  fasta=2')
# args.add_argument('--cleanfile', default='cleandf file path')
args.add_argument('--reward', default=0.1, help='reward for new reactions to update predscore')
args.add_argument('--iter', default = 3, help='The number of iterations to run the funcarve algorithm')

args.add_argument('--block_flage', default = 1, help='decrease block reactions: 1, ignore block reactions: 0')
args.add_argument('--flux_flage', default = 0, help='reward changed by flux of the reactions: 1, ignore reactions fluxes: 0')
args.add_argument('--media', default='rich', help='List of metabolites composing the media condition. Not required.')
args.add_argument('--tasks', default=[], help='List of metabolic tasks. Not required.')
args.add_argument('--org', default='default', help='KEGG organism code. Not required.')
args.add_argument('--min_frac', default=0.01, help='Minimum objective fraction required during gapfilling')
args.add_argument('--max_frac', default=0.5, help='Maximum objective fraction allowed during gapfilling')

args.add_argument('--threshold', default = 8, help='The cutoff value for the EC prediction score')
args.add_argument('--upper', default=15, help='Upper threshold for predscore')
args.add_argument('--lower', default=5, help='Lower threshold for predscore')
args.add_argument('--maxweight', default=100, help='Maximum weight for reactions')
args.add_argument('--minweight', default=0.0, help='Minimum weight for reactions')

# args.add_argument('--threshold', default = 0.5, help='The cutoff value for the EC prediction score')
# args.add_argument('--upper', default=15, help='Upper threshold for predscore')
# args.add_argument('--lower', default=5, help='Lower threshold for predscore')
# args.add_argument('--maxweight', default=100, help='Maximum weight for reactions')
# args.add_argument('--minweight', default=0.0, help='Minimum weight for reactions')

args.add_argument('--gram', default='none', help='Type of Gram classificiation (positive or negative)')
args.add_argument('--out', default='default', help='Name of output GENRE file')
args.add_argument('--name', default='default', help='ID of output GENRE')
args.add_argument('--cpu', default=1, help='Number of processors to use')
args.add_argument('--gapfill', default='yes', help='gapfill your model?')
args.add_argument('--exchange', default = 1, help='open exchange: 1, shut down exchange: 0')
args.add_argument('--test', default = 'no', help='do you want to perform the test suite?')

args = args.parse_args()


with open('../data/biggr2ec.pkl', 'rb') as f:
    biggr2ec = pickle.load(f)
with open('biggec2r.pkl', 'rb') as f:
    biggec2r = pickle.load(f)
print("biggr2ec dictionary loaded successfully.")
with open('../data/seedr2ec.pkl', 'rb') as f:
    seedr2ec = pickle.load(f)
with open('../data/seedec2r.pkl', 'rb') as f:
    seedec2r = pickle.load(f)
print("seedr2ec dictionary loaded successfully.")

  
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
   
#----------------------------------------------------------------------------------------------------------------------#
if __name__ == "__main__":
    # Process input settings
    input_file = str(args.input_file)
    # cleanfile = str(args.cleanfile)
    out_file = str(args.out)
    file_type = int(args.file_type)
    name = str(args.name)
    # org = str(args.org)
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
    
    # upper=15,lower=5,maxweight=100,minweight=0.0
    upper = float(args.upper)
    lower = float(args.lower)
    maxweight = float(args.maxweight)
    minweight = float(args.minweight)

    ## print all input arguments in condense form
    print('all input arguments:\n',args,sep='')


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

    # if org != 'default':
    #     print('Including additional genes from KEGG genome of', org)

    # Maximum fraction should not be too high, otherwise the gapfiller adds too many reactions
    print('Using minimum objective flux fraction of', min_frac,'and maximum fraction of', max_frac)

    if processors > cpu_count():
        print('WARNING: Requested more processors than are available. Using maximum of', cpu_count())
        processors = cpu_count()
    print('Using', processors, 'processor(s)\n')


    # Load databases
    with open('../data/universal.pickle', 'rb') as f:
        universal = pickle.load(f)
    with open('../data/gene_names.pickle', 'rb') as f:
        gene_names = pickle.load(f)

    for i in range(1,int(inter)+1):
        print('*'*50)
        print('Inter:',i)
        print('inter=',inter)
        print('*'*50)
        if i == 1:
            if file_type ==2:
                fasta_file = input_file
                clean_file = run_clean(fasta_file)
            pr2ec,ec2pr,predscore,allec2pr = read_cleandf_withscore(clean_file,threshold=threshold)

            ## save the original predscore
            with open(f'../data/tmp/{name}_t{threshold}_ORIpredscore_{i-1}.pkl', 'wb') as f:
                pickle.dump(predscore, f)
            print(f'saved original predscore from CLEAN prediction',flush=True)

        elif i > 1:
            if len(updateprs) == 0:
                print('No updateprs')
                print(f'End No more iteration no more {i}############')
                break
            predscore = newpredscore

        updateprs = []
        universal_scoredict,rxns = clean2seedr(predscore,threshold=threshold)
        print('Draft clean_to_rxns has', len(rxns), 'reactions',flush=True)
        newpredscore = deepcopy(predscore)
        draft_genre = create_model(rxns, universal, new_id)
        print('Draft draft_genre has', len(draft_genre.reactions), 'reactions',flush=True)
        print('Draft draft_genre has', len(draft_genre.metabolites), 'metabolites',flush=True)
        print('Draft draft_genre has', len(draft_genre.genes), 'genes',flush=True)

        if block_flage == 1:
            block_rxns = cobra.flux_analysis.find_blocked_reactions(draft_genre)
            print('Draft draft_genre has', len(block_rxns),'block reactions',flush=True)
            newpredscore, updateprs, allec2pr = update_predscore_block(block_rxns,allec2pr,newpredscore,updateprs,reward,threshold)

        draft_genre = add_names(draft_genre, gene_names)
        
        
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
            # simulated_media = [
            media = [
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
        print('Using media condition:', len(media),'metabolites for exchange',flush=True)
        
        # Set media condition
        if len(media) != 0:
            media_condition = set(['EX_' + cpd for cpd in media])
            universal_reactions = set([x.id for x in universal.reactions])
            for rxn in universal_reactions:
                if rxn.startswith('EX_') == True:
                    universal.reactions.get_by_id(rxn).bounds = (0, 1000.0)
                if rxn in media_condition:
                    universal.reactions.get_by_id(rxn).bounds = (-1000.0, 10000)

        draft_reactions = set([x.id for x in draft_genre.reactions])
        draft_metabolites = set([x.id for x in draft_genre.metabolites])

        # Gapfill new model
        if gapfill == 'yes':
            ##############################################
            print('>>Identifying new metabolism (Step 1 of 2)...')    
            ##############################################
            new_reactions ,coefficientDict,new_rxn_flux = weighted_find_reactions(draft_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 1, file_type,universal_scoredict,upper,lower,maxweight,minweight)
            print('STEP1:Gapfilling new reaction:',len(new_reactions),flush=True)
            newpredscore, updateprs, allec2pr = update_predscore_flux(new_rxn_flux,allec2pr,newpredscore,updateprs,reward,threshold,flux_flage)
            
            ## remove unsatisfactory reactions
            if gram_type == 'positive':
                new_reactions = set([x for x in new_reactions if x not in ['GmNeg_cellwall','biomass_GmNeg']])
            elif gram_type == 'negative':
                new_reactions = set([x for x in new_reactions if x not in ['GmPos_cellwall','biomass_GmPos']])
            print('STEP1:Gapfilling new reactions af remove gram dismatch reaction',len(new_reactions),flush=True)
            filled_genre = gapfill_model(draft_genre, universal, new_reactions, universal_obj, 1)


            ##############################################
            print('>>Identifying new metabolism (Step 2 of 2)...')
            ##############################################
            filled_genre = set_base_inputs(filled_genre, universal)
            media_reactions,coefficientDict,media_flux = weighted_find_reactions(filled_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 2, file_type,universal_scoredict,upper,lower,maxweight,minweight)
            print('STEP2:Gapfilling new media reaction:',len(media_reactions),flush=True)
            newpredscore,updateprs,allec2pr = update_predscore_flux(media_flux,allec2pr,newpredscore,updateprs,reward,threshold,flux_flage)
            
            if gram_type == 'positive':
                media_reactions = set([x for x in media_reactions if x not in ['GmNeg_cellwall','biomass_GmNeg']])
            elif gram_type == 'negative':
                media_reactions = set([x for x in media_reactions if x not in ['GmPos_cellwall','biomass_GmPos']])
            print('STEP2:Gapfilling new reactions af remove gram dismatch reaction',len(media_reactions),flush=True) 

            final_genre = gapfill_model(filled_genre, universal, media_reactions, universal_obj, 2)
            
        else:
            final_genre = draft_genre
            
        final_genre = add_annotation(final_genre, gram_type)
        print('# Correct exchanges and check new model')
        if exchange_arg == 0:
            for exch in final_genre.exchanges: exch.bounds = (0., 0.)
        else:
            for exch in final_genre.exchanges: exch.bounds = (-1000., 1000.)
        for rxn in final_genre.reactions:
            if 'Exchange reaction for' in rxn.name:
                rxn.name = list(rxn.metabolites)[0].name + ' exchange'
        biomass = checkModel(draft_reactions, draft_metabolites, final_genre)
        
        print('>>> Write model to sbml')
        input_file = input_file.split('/')[-1] # write to working directory
        
        if file_type == 1:
            if new_id != 'default':
                out_file = input_file.rstrip('maxsep_df.pkl') + new_id + f'I{i}'+ '.sbml'
            else:
                out_file = input_file.rstrip('maxsep_df.pkl') + 'enzbuild' +  f'I{i}'+'.sbml'
        elif file_type == 2:
            if new_id != 'default':
                out_file = input_file.rstrip('.fasta') + new_id + f'I{i}'+ '.sbml'
            else:
                out_file = input_file.rstrip('.fasta') + 'enzbuild' +  f'I{i}'+'.sbml'
            out_file = '../data/result/' + out_file
        print('\n>>Saving new GEM to', out_file, '\n')

        ## save the dictionary
        thr = str(threshold).split('.')[-1]
        if float(threshold) >= 1.0:
            thr = str(float(threshold)).split('.')[0]
        else:
            thr = str(threshold).split('.')[-1]



        with open(f'../data/tmp/{name}_t{thr}_differ_predscore_{i}.pkl', 'wb') as f:
            pickle.dump(newpredscore, f)
        print(f'saved differ predscore {i}')
        with open(f'../data/tmp/{name}_t{thr}_updateprs_{i}.pkl', 'wb') as f:
            pickle.dump(updateprs, f)
        print(f'saved updateprs {i}')

        cobra.io.write_sbml_model(final_genre, out_file)
        print('*'*50)
        print('End of Inter:',i)
        print('*'*50)
