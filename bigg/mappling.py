import pandas as pd
import json
from utils import *

metaf = '/ibex/user/niuk0a/funcarve/cobra/bigg/bigg_models_metabolites.txt'
reactionf = '/ibex/user/niuk0a/funcarve/cobra/bigg/bigg_models_reactions.txt'
modelf = '/ibex/user/niuk0a/funcarve/cobra/bigg/universal_model.json'
r2gf = '/ibex/user/niuk0a/funcarve/cobra/bigg/rhea2go'
ec2gf = '/ibex/user/niuk0a/funcarve/cobra/bigg/ec2go'
kr2gf = '/ibex/user/niuk0a/funcarve/cobra/bigg/kegg_reaction2go'
# /ibex/user/niuk0a/funcarve/cobra/bigg/kegg_reaction2go
# rhea2go = read_rhea2go(r2gf)
## reactions 
reactions = pd.read_csv(reactionf, sep='\t')

# print(reactions.columns)
# print(list(reactions.name.values)[:10])
# print(list(reactions.database_links.values)[:2])
# for i in list(reactions.database_links.values)[:2]:
#     print('reaction:',i)



# RHEA 2 GO
# f = '/ibex/user/niuk0a/funcarve/cobra/bigg/rhea2go'

count = 0
for link in list(reactions.database_links.values):
    link = database_links_reformat(link)
    # print('link',link)
    ids, rhea_ids, mnxs, seeds, biocyc, ec, keggr = links_to_id(link)
    # print('ids:',ids)
    # print('rhea_ids:',rhea_ids)
    # print('mnxs:',mnxs)
    # print('seeds:',seeds)
    # print('biocyc:',biocyc)
    # rhea2go = read_rhea2go(r2gf)
    ec2go = read_ec2go(ec2gf)
    # keggr2go = read_keggr2go(kr2gf)

    # go = get_rhea_go(rhea2go, rhea_ids)
    # go_kegg = get_keggr_go(keggr2go, keggr)
    go_ec = get_ec_go(ec2go, ec)
    ## go = get_mnx_go(mnxs)
    ## go = get_seed_go(seeds)
    ## go = get_biocyc_go(biocyc)


    # combine go and go_ec
    # print('len(go)', len(go))
    # print('len(go_ec)', len(go_ec))
    # go = go.union(go_ec)
    # go = go.union(go_kegg)
    go = go_ec # only use ec2go
    # print('only use ec2go')
    # print('len(go)', len(go))
    # print('-'*10)
    if go != set():
        count += 1
    # matching_reactions = reactions[reactions.database_links.apply(lambda x: database_links_reformat(x) == link)]
        # print('-'*10)
        # print('r:', 'maps to', go)
        # print('-'*10) 
print('only use ec2go')

print(f'count:{count}, no go:{len(reactions)-count}, total:{len(reactions)}')
# Rhea count:1972, no go:26329, total:28301
# which means that 26329 reactions are not mapped to GO from RHEA

# then RHEA+EC count:3897, no go:24404, total:28301
# then RHEA+EC+KEGG count:3899, no go:24402, total:28301
#only use ec2go count:3566, no go:24735, total:28301





# MetaNetX 2 GO
# SEED 2 GO



# metabolites = pd.read_csv(metaf, sep='\t')
# reactions = pd.read_csv(reactionf, sep='\t')

# print(metabolites.columns)
# print(reactions.columns)

# print(metabolites.head())
# print(reactions.head())