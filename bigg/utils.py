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

