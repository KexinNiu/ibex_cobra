import cobra
import pickle

from utils import read_fasta
import pandas as pd
import numpy as np
# def writefasta(names,seqs,filename):
#     with open(filename, 'w') as f:
#         for i in range(len(names)):
#             name = names[i]
#             if len(name.split('[gene=')) >1:
#                 name = name.split('[gene=')[1].split(']')[0]
#             else:
#                 name = name.split('[locus_tag=')[1].split(']')[0]
#             f.write('>'+name+'\n')
#             f.write(seqs[i]+'\n')
# names,seqs = read_fasta('/ibex/user/niuk0a/funcarve/cobra/NC_000962.3_sequence_inj661.txt')
# writefasta(names,seqs,'/ibex/user/niuk0a/funcarve/cobra/NC_000962.3_sequence_inj661.fasta')
# names,seqs = read_fasta('/ibex/user/niuk0a/funcarve/cobra/CP000647.1_sequence_iyj1228.txt')
# writefasta(names,seqs,'/ibex/user/niuk0a/funcarve/cobra/CP000647.1_sequence_iyj1228.fasta')


# exit()

with open('/ibex/user/niuk0a/funcarve/cobra/Unique_ModelSEED_Reaction_ECs.txt','r') as f:
    r2ec={}
    ec2r={}
    for line in f:
        if line.startswith('ModelSEED'):
            continue
        
        else:
            item = line.strip('\n')
            items = item.split('\t')
            # r2ec[items[0]] = items[1]
            # ec2r[items[1]] = items[0]
            try:
                r2ec[items[0]].append(items[1])
            except KeyError:
                r2ec[items[0]] = [items[1]]
            try:    
                ec2r[items[1]].append(items[0])
            except KeyError:
                ec2r[items[1]] = [items[0]]
    #save the dictionary
    with open('seedr2ec.pkl', 'wb') as f:
        pickle.dump(r2ec, f)
    with open('seedec2r.pkl', 'wb') as f:
        pickle.dump(ec2r, f) 
        
exit()
f = '/ibex/user/niuk0a/funcarve/cobra/SEED2VMH_translation.csv'
ms2bigg={}
bigg2ms={}
with open(f, 'r') as file:
    for line in file:
        line = line.strip('\n').split(',')
        ms2bigg[line[0]] = line[1]
        

universal = cobra.io.load_json_model("/ibex/user/niuk0a/funcarve/cobra/bigg/universal_model_cobrapy.json")
print('len(universal.reactions)',len(universal.reactions))

## remove reactions
f = '/ibex/user/niuk0a/funcarve/cobra/recon_remove.txt'
remove_bigg = []
with open(f, 'r') as file:
    next(file)
    for line in file:
        line = line.strip('\n')
        id = line.split('_')[0]
        try:
            biggid = ms2bigg[id]
            remove_bigg.append(universal.reactions.get_by_id(biggid))
        except:
            # print('id:',id)
            continue
print('len(remove_bigg)',len(remove_bigg))
        
    
        

# universal = cobra.io.read_sbml_model("/ibex/user/niuk0a/funcarve/cobra/bigg/bigg_universe.xml")

universal.remove_reactions(remove_bigg)
print('len(universal.reactions)',len(universal.reactions))
universal.repair()
#save universal
cobra.io.write_sbml_model(universal, 'universal_model_cobrapy_filter.xml')
# exit()

# filename='/ibex/user/niuk0a/funcarve/cobra/universal_recon.pickle'
# with open(filename, 'rb') as f: 
#     reconuniversal = pickle.load(f)

# # filter universal
# print('len(reconuniversal.reactions)',len(reconuniversal.reactions))
# uuu  = len(universal.reactions)
# newrecon = cobra.Model('newrecon')
# urm = []
# for r in reconuniversal.reactions:
#     msid = r.id
#     # print('msid:',msid)
#     if msid.startswith('SNK_'):
#         msid = msid[4:-2]
#         msid = 'EX_'+msid+'(e)'
#     elif msid.startswith('rxn'):
#         msid = msid[:-2]  
#     try:
#         biggid = ms2bigg[msid]
#         print('msid:',msid,'biggid:',biggid)
#         ur = universal.reactions.get_by_id(biggid)
#         ur.id = biggid
#     except:
#         urm.append(r)
# print('len(urm)',len(urm))
# universal.remove_reactions(urm)
# universal.repair()
# #save universal
# print('len(urm)',len(urm))
# cobra.io.write_sbml_model(universal, 'filter_bigg_universe.xml')
# print('len(reconuniversal.reactions)',len(reconuniversal.reactions))

# print('len(universal.reactions)','before',uuu,'now',len(universal.reactions))  


####model 2 modelseed 
# f = '/ibex/user/niuk0a/funcarve/cobra/iAF987_555_0.05_rm.xml'
# model = cobra.io.read_sbml_model(f)
# for r in model.reactions:
#     if r.id in ms2bigg:
#         r.id = ms2bigg[r.id]
