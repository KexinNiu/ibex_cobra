import requests
import pandas as pd

# Base URL for the BiGG Models API
base_url = 'http://bigg.ucsd.edu/api/v2'
# Get a list of models
# curl 'http://bigg.ucsd.edu/api/v2/models'
# Get a list of model reactions
# curl 'http://bigg.ucsd.edu/api/v2/models/iND750/reactions'
# curl 'http://bigg.ucsd.edu/api/v2/models/iMM904/genes'
# curl 'http://bigg.ucsd.edu/api/v2/models/iMM904/genes/Q0045'

def get_all_models():
    url = f'{base_url}/models'
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()['results']
    else:
        print(f"Error fetching models: {response.status_code}")
        return []

# save list to file
def make_df_bigg_curlresult(modellist):
    bigg_id = []
    gene_count = []
    reaction_count = []
    organism = []
    metabolite_count = []

    for item in modellist:
        item = dict(item)
        bigg_id.append(item['bigg_id'])
        organism.append(item['organism'])
        gene_count.append(item['gene_count'])
        reaction_count.append(item['reaction_count'])
        metabolite_count.append(item['metabolite_count'])
    df = pd.DataFrame({'bigg_id': bigg_id, 'organism': organism, 'gene_count': gene_count, 'reaction_count': reaction_count, 'metabolite_count': metabolite_count})
    # df.to_csv("bigg_models.csv", index=False)
    return df

# get all models
# modellist = get_all_models()
# df = make_df_bigg_curlresult(modellist)
# print(df.head())
# #save the df to a csv file
# df.to_csv("bigg_models.csv", index=False)

################################################################################
# get genes for each model
# curl 'http://bigg.ucsd.edu/api/v2/models/e_coli_core'

def get_gene_list_for_model(model_id):
    url = f'{base_url}/models/{model_id}/genes'
    response = requests.get(url)
    if response.status_code == 200:
        gene_results = response.json()['results']
        gene_list_biggid = [item['bigg_id'] for item in gene_results]
        gene_list_name = [item['name'] for item in gene_results]

        return gene_list_biggid, gene_list_name
    else:
        print(f"Error fetching genes: {response.status_code}")
        return [], []

def get_gene_meta(model_id):
    gene_list_biggid, _ = get_gene_list_for_model(model_id)
    bigg_id = []
    name = []
    strand = []
    leftpos = []
    genome_ref_string = []
    database_links = []
    rightpos = []
    protein_sequence = []
    model_bigg_id = []
    reactions = []
    old_identifiers = []
    mapped_to_genbank = []
    chromosome_ncbi_accession = []

    for gene_biggid in gene_list_biggid:
        url = f'{base_url}/models/{model_id}/genes/{gene_biggid}'
        response = requests.get(url)
        if response.status_code == 200:
            gene_meta = response.json()
            bigg_id.append(gene_meta['bigg_id'])
            name.append(gene_meta['name'])
            strand.append(gene_meta['strand'])
            leftpos.append(gene_meta['leftpos'])
            genome_ref_string.append(gene_meta['genome_ref_string'])
            database_links.append(gene_meta['database_links'])
            rightpos.append(gene_meta['rightpos'])
            protein_sequence.append(gene_meta['protein_sequence'])
            model_bigg_id.append(gene_meta['model_bigg_id'])
            reactions.append(gene_meta['reactions'])
            old_identifiers.append(gene_meta['old_identifiers'])
            mapped_to_genbank.append(gene_meta['mapped_to_genbank'])
            chromosome_ncbi_accession.append(gene_meta['chromosome_ncbi_accession'])
        else:
            print(f"Error fetching gene meta: {response.status_code}")

    df = pd.DataFrame({'bigg_id': bigg_id, 
                        'name': name, 
                        'strand': strand, 
                        'leftpos': leftpos, 
                        'genome_ref_string': genome_ref_string, 
                        'database_links': database_links, 
                        'rightpos': rightpos, 
                        'protein_sequence': protein_sequence, 
                        'model_bigg_id': model_bigg_id, 
                        'reactions': reactions, 
                        'old_identifiers': old_identifiers, 
                        'mapped_to_genbank': mapped_to_genbank, 
                        'chromosome_ncbi_accession': chromosome_ncbi_accession})
    print(df.head(),flush=True)
    return df

modeldf = pd.read_csv("bigg_models.csv")
print('done reading modeldf',flush=True)

proteinname = []
proteinsequence = []
# if 

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--part', type=int, help='Part of the models to process')

args = parser.parse_args()


interval = 108
part = args.part
start = int(interval / 4) * part
end = int(interval / 4) * (part+1)
if part == 3:
    end = interval+1
fastafile = open(f"bigg_genes_part_{part}.fasta", "w")
print(f"#########################\nstart {part}\nstart {start}\tend {end}\n#########################\n",flush=True)

finalmetadf = pd.DataFrame()
metadf = pd.DataFrame()

for index, row in modeldf.iterrows():
    if index < start or index >= end:
        continue
    if index >= end:
        break

    print(f"start {index}",flush=True)
    model_id = row['bigg_id']
    print(f"start {model_id}",flush=True)
    df = get_gene_meta(model_id)
    df.to_pickle(f"bigg_genes_meta_{model_id}.pkl")
    metadf = pd.concat([metadf, df])
    metadf.to_pickle(f"bigg_genes_meta_tmp_{part}.pkl")

    print('done getting gene meta',part,flush=True)
    for i, row in df.iterrows():
        if row['name'] == 'None':
            print(f"skip {row['model_bigg_id']} None")
            continue
        tmpname = str(row['model_bigg_id'])+'_'+str(row['name']) 
        if row['protein_sequence'] in proteinsequence or row['protein_sequence'] == '':
            if row['protein_sequence'] == '':
                print(f"skip {tmpname} empty")
            else:
                print(f"skip {tmpname} already in")
            continue
        else:
            proteinname.append(tmpname) 
            proteinsequence.append(row['protein_sequence'])
            fastafile.write(f">{tmpname}\n{row['protein_sequence']}\n")
    print(f"done {model_id}",'{}/{}'.format(index+1, len(modeldf)),flush=True)
metadf.to_pickle(f"bigg_genes_meta_final_{part}.pkl")
fastafile.close()
print('done writing fasta',flush=True)
metadf.to_pickle(f"bigg_genes_meta_{part}.pkl")
print('done writing meta',flush=True)

# f = '/ibex/user/niuk0a/funcarve/cobra/bigg/bigg_genes_meta_tmp_0.pkl'
f0='/ibex/user/niuk0a/funcarve/cobra/bigg/bigg_genes_part_0.fasta'
f1='/ibex/user/niuk0a/funcarve/cobra/bigg/bigg_genes_part_1.fasta'
f2='/ibex/user/niuk0a/funcarve/cobra/bigg/bigg_genes_part_2.fasta'
f3='/ibex/user/niuk0a/funcarve/cobra/bigg/bigg_genes_part_3.fasta'
n0,s0 = read_fasta(f0)
n1,s1 = read_fasta(f1)
n2,s2 = read_fasta(f2)
n3,s3 = read_fasta(f3)

n0,s0 = read_fasta(f0)
n1,s1 = read_fasta(f1)
n2,s2 = read_fasta(f2)
n3,s3 = read_fasta(f3)
alls = set(s0)
alls.update(s1)
alls.update(s2)
alls.update(s3)
len(alls)
41110

def dic_d(n,s,d):
    if d == {}:
        for i in range(0,len(n)):
            d[s[i]] = [n[i]]
    else:
        for i in range(0,len(n)):
            try:
                d[s[i]].append(n[i])
            except:
                d[s[i]] = [n[i]]
    return d

d = dic_d(n0,s0,{})
d1= dic_d(n1,s1,d)
d2= dic_d(n2,s2,d1)
d3= dic_d(n3,s3,d2)
len(d3)
41113
d3.pop('None')
['e_coli_core_None', 'iEC1372_W3110_S0001', 'iECSP_1301_None', 'iNJ661_None']
41112

with open('bigg_sequences.pickle', 'wb') as handle:
    pickle.dump(d3, handle, protocol=pickle.HIGHEST_PROTOCOL)


