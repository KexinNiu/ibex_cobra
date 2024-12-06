refs = '/ibex/user/niuk0a/funcarve/cobra/NC_000913.3.txt'
cleanf = '/ibex/user/niuk0a/funcarve/cobra/NC_000913.3_ecoli_maxsep.csv'
unif='/ibex/user/niuk0a/funcarve/cobra/NC_000913.3_uniprotkb_taxonomy_id_83333_2024_10_23.tsv'



refs = '/ibex/user/niuk0a/funcarve/cobra/aaseq_CP000148.1.txt'
cleanf = '/ibex/user/niuk0a/funcarve/cobra/aaseq_CP000148.1_maxsep.csv'
unif='/ibex/user/niuk0a/funcarve/cobra/uniprotkb_taxonomy_id_269799_AND_review_2024_10_27.tsv'
import pandas as pd
def read_refseq_fasta(fasta_file):
    seq=''
    names =[]
    seqs = []
    genes = []
    db_xref = []
    protein_id = []
    locustag = []
    with open(fasta_file, 'r') as inFile:
        for line in inFile:
            if line.startswith('>'):
                name = line.strip('\n').split('>')[1]
                cleanname = name.split('] [')[0].split('lcl|')[1]
                names.append(cleanname)
                items = name.split(' ')
                for item in items:
                    if item.startswith('['):
                        item = item[1:-1]
                        if '=' in item:
                            key,val = item.split('=')
                        else:
                            continue
                        if key == 'gene':
                            genes.append(val)
                        elif key == 'locus_tag':
                            locustag.append(val)
                        elif key =='db_xref':
                            if val.startswith('UniProtKB/Swiss-Prot:'):
                                val = val.split(':')[-1]
                            else:
                                # print('val:',val)
                                val = val.split(':')[-1]
                            db_xref.append(val)
                        elif key =='protein_id':
                            protein_id.append(val)
                if len(names) == len(genes) == len(db_xref) == len(protein_id):
                    pass
                else:
                    ## add a none val to the list
                    if len(names) > len(genes):
                        genes.append('None')
                    if len(names) > len(db_xref):
                        db_xref.append('None')
                    if len(names) > len(protein_id):
                        protein_id.append('None')           
                # try:
                #     name = line.strip('\n').split('[gene=')[1].split(']')[0]
                # except IndexError:
                #     name = line.strip('\n').split('[locus_tag=')[1].split(']')[0]   
                
                
                if seq == '':
                    continue
                else:
                    seqs.append(seq)
                    seq = ''    
            else:
                seq = seq + line.strip('\n')
        seqs.append(seq)

    return names,seqs,genes,db_xref,protein_id,locustag

def read_clean_withscore(input_file,threshold=0.8):
    # print('threrhold-->',threshold)
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
                    if dis >= -0.0001:
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

def read_unif(f):
    df = pd.read_csv(f,sep='\t')
    # only take where Reviewed == reviewed
    r_df = df[df['Reviewed']=='reviewed'] 
    #  take columns ['Entry','Protein names','Gene names','Organism','EC number','Status']
    r_df = r_df[['Entry','Gene Names','EC number']]
    return r_df

def check_inter(uniprotecs,cleanecs,counts,total):
    # re = [0,0,0,0]
    # ptoto =[0,0,0,0]
    flage=False
    if ';' in uniprotecs:
        uniprotecs = uniprotecs.split('; ')
    for e in uniprotecs:
        level = e.count('-')
        for ec in cleanecs:
            total[level] +=1
            if ec.startswith(e.split('-')[0]):
                counts[level] +=1
                flage=True     
    return counts,total,flage
            
def print_sta(r_df,db_xref,predscore,thres):   
    count =[0,0,0,0]
    total = [0,0,0,0]
    totalcount = 0
    number = 0
    # thres =0.6
    for p in db_xref:
        cleanname = names[db_xref.index(p)]
        try:
            cleanpred = predscore[cleanname]
        except:
            # print('skip protein:',p)
            continue
        
        unip = r_df[r_df['Entry']==p]
        if unip.empty:
            continue
        if unip['EC number'].isnull().values.any():
            continue
        
        # print('protein:',p)
        # print(cleanpred,unip)
        uniecs = unip['EC number'].values.tolist()[0]
        # cleanecs = set([i.split(':')[1] for i in cleanpred.keys() if cleanpred[i] >= thres])
        cleanecs = set([i for i in cleanpred.keys() if cleanpred[i] >= thres])
        
        number+=1
        count,total,flage = check_inter(uniecs,cleanecs,count,total)
        
        if flage:
            totalcount+=1
        
    print('threshold=',thres)  
    # print('digits 1:',count[0],total[0],count[0]/total[0])
    # print('digits 2:',count[1],total[1],count[1]/total[1])
    # print('digits 3:',count[2],total[2],count[2]/total[2])
    # print('digits 4:',count[3],total[3],count[3]/total[3])
    # print('total:',totalcount,len(db_xref),totalcount/len(db_xref))

    print('digits 1:',count[0],total[0])
    print('digits 2:',count[1],total[1])
    print('digits 3:',count[2],total[2])
    print('digits 4:',count[3],total[3])
    print('total:',totalcount,len(db_xref),totalcount/len(db_xref))

    # print('level 1',count[0]/total[0],sep='|')
    # print('level 2',count[1]/total[1],sep='|')
    # print('level 3',count[2]/total[2],sep='|')
    # print('level 4',count[3]/total[3],sep='|')
    # print('total',totalcount/number,sep='|')
    return 



names, _, genes, db_xref, protein_id,locustag = read_refseq_fasta(refs)
print('len of all output:',len(names),len(genes),len(db_xref),len(protein_id),len(locustag))
pr2ec,predscore = read_clean_withscore(cleanf)
print(names[:10])
print('genes:',genes[:10])
print('protein_id:',protein_id[:10])
cout = [i for i in db_xref if i != 'None']
print('protein number:',len(cout),cout)
print('locustag:',locustag[:10])
print('db_xref:',len(db_xref),db_xref[:10])
print('pr2ec:',len(pr2ec))

# print('predscore:',predscore)
# newprediscoref='/ibex/user/niuk0a/funcarve/cobra/ecoliINTER_newpredscore_8.pkl'
newprediscoref='/ibex/user/niuk0a/funcarve/cobra/iaf987INTER_newpredscore_3.pkl'

newprediscore = pd.read_pickle(newprediscoref)
# print('newprediscore:',newprediscore)
# r_df = read_unif(unif)
# for th in [0.4,0.5,0.6,0.7,0.8,0.9]:
#     print_sta(r_df,db_xref,newprediscore,th)

# print('original predscore:_________________________')
# for th in [0.4,0.5,0.6,0.7,0.8,0.9]:
#     print_sta(r_df,db_xref,predscore,th)



