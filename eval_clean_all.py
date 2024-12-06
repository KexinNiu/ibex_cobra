from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.metrics import precision_score, recall_score, \
    roc_auc_score, accuracy_score, f1_score, average_precision_score
import pandas as pd
import csv
import numpy as np
import os
from Bio import SeqIO


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

def get_pred_labels(out_filename, pred_type="_maxsep"):
    file_name = out_filename+pred_type
    result = open(file_name+'.csv', 'r')
    csvreader = csv.reader(result, delimiter=',')
    pred_label = []
    pred_names = []
    pred_all_label=set()
    for row in csvreader:
        preds_ec_lst = []
        preds_with_dist = row[1:]
        pred_names.append(row[0])
        for pred_ec_dist in preds_with_dist:
            # get EC number 3.5.2.6 from EC:3.5.2.6/10.8359
            ec_i = pred_ec_dist.split(":")[1].split("/")[0]
            preds_ec_lst.append(ec_i)
            pred_all_label.add(ec_i)
        pred_label.append(preds_ec_lst)
    # print('>get_pred_labels:',len(pred_label),len(pred_names),len(pred_all_label))
    return pred_label, pred_names,pred_all_label

def get_pred_probs(out_filename, pred_type="_maxsep"):
    file_name = out_filename+pred_type
    result = open(file_name+'.csv', 'r')
    csvreader = csv.reader(result, delimiter=',')
    pred_probs = []
    for row in csvreader:
        preds_ec_lst = []
        preds_with_dist = row[1:]
        # probs = torch.zeros(len(preds_with_dist))
        probs = np.zeros(len(preds_with_dist))
        count = 0
        for pred_ec_dist in preds_with_dist:
            # get EC number 3.5.2.6 from EC:3.5.2.6/10.8359
            ec_i = float(pred_ec_dist.split(":")[1].split("/")[1])
            probs[count] = ec_i
            #preds_ec_lst.append(probs)
            count += 1
        # sigmoid of the negative distances 
        # probs = (1 - torch.exp(-1/probs)) / (1 + torch.exp(-1/probs))
        probs = (1 - np.exp(-1/probs)) / (1 + np.exp(-1/probs))
        # probs = probs/torch.sum(probs)
        probs = probs/np.sum(probs)
        pred_probs.append(probs)
    return pred_probs

def get_eval_metrics(pred_label, pred_probs, true_label, all_label):
    mlb = MultiLabelBinarizer()
    mlb.fit([list(all_label)])
    n_test = len(pred_label)
    pred_m = np.zeros((n_test, len(mlb.classes_)))
    true_m = np.zeros((n_test, len(mlb.classes_)))
    # for including probability
    pred_m_auc = np.zeros((n_test, len(mlb.classes_)))
    label_pos_dict = get_ec_pos_dict(mlb, true_label, pred_label)
    for i in range(n_test):
        pred_m[i] = mlb.transform([pred_label[i]])
        true_m[i] = mlb.transform([true_label[i]])
         # fill in probabilities for prediction
        labels, probs = pred_label[i], pred_probs[i]
        for label, prob in zip(labels, probs):
            if label in all_label:
                pos = label_pos_dict[label]
                pred_m_auc[i, pos] = prob
    pre = precision_score(true_m, pred_m, average='weighted', zero_division=0)
    rec = recall_score(true_m, pred_m, average='weighted')
    f1 = f1_score(true_m, pred_m, average='weighted')
    roc = roc_auc_score(true_m, pred_m_auc, average='weighted')
    acc = accuracy_score(true_m, pred_m)
    return pre, rec, f1, roc, acc

def get_pred_labels_prc(out_filename, cutoff, pred_type="_maxsep"):
    file_name = out_filename+pred_type
    result = open(file_name+'.csv', 'r')
    csvreader = csv.reader(result, delimiter=',')
    pred_label = []
    pred_all_label=set()
    for row in csvreader:
        preds_ec_lst = []
        preds_with_dist = row[1:]
        for pred_ec_dist in preds_with_dist:
            # get EC number 3.5.2.6 from EC:3.5.2.6/10.8359
            ec_i = pred_ec_dist.split(":")[1].split("/")[0]
            if int(pred_ec_dist.split(":")[1].split("/")[1]) <= cutoff:
                preds_ec_lst.append(ec_i)
                pred_all_label.add(ec_i)
        pred_label.append(preds_ec_lst)
    return pred_label ,pred_all_label

def get_true_labels(file_name,pref,pred_names,pred_label,pred_probs):
    
    gbkf = pref+'.gb'
    r_df,df,lt2ec,lt2oldlt, oldlt2lt = read_unif_genebk(file_name,gbkf)
    all_label = set()
    true_label_dict = {}
    count=0
    rm_index = []
    goinfo={}
    true_label = []
    
    currpred_names = []
    currpred_label = []
    currpred_probs = []
    if not lt2ec == {}:
    #########################################################################
    #################### WITH GENEBANK + UNIPROT ANNOTATION #################
        for i in range(0,len(pred_names)):
            name = pred_names[i]
            try:
                oldname = lt2oldlt[name]
            except:
                continue
            remove = True
            for on in oldname:
                unip = r_df[r_df['Gene Names'].str.contains(on, na=False)]
                if len(unip) == 0:
                    continue
                elif unip['EC number'].isnull().values.any():
                    continue
                    goinfo[name] = unip['Gene Ontology (GO)']
                elif not unip['EC number'].isnull().values.any():
                    remove = False
                    uniecs = unip['EC number'].values.tolist()[0]
                    true_ec_lst = uniecs.split('; ')
                    for ec in true_ec_lst:
                        all_label.add(ec)
                    true_label_dict[name] = true_ec_lst
                    break
                elif len(lt2ec[on]) > 0:
                    remove = False
                    true_ec_lst = lt2ec[on]
                    for ec in true_ec_lst:
                        all_label.add(ec)
                    true_label_dict[name] = true_ec_lst
            if remove:
                rm_index.append(i)
            else:
                currpred_names.append(name)
                currpred_label.append(pred_label[i])
                currpred_probs.append(pred_probs[i])
                true_label.append(true_label_dict[name])
    #########################################################################

    ###################################################################
    #################### WITH ONLY UNIPROT ANNOTATION #################
    else:
        for i in range(0,len(pred_names)):
            name = pred_names[i]
            unip = r_df[r_df['Gene Names'].str.contains(name, na=False)]
            if len(unip) == 0:
                # print('no uniprot:',name,uniname)
                    rm_index.append(i)
                    continue
            elif unip['EC number'].isnull().values.any():
                rm_index.append(i)
                goinfo[name] = unip['Gene Ontology (GO)']
                continue
            else:
                uniecs = unip['EC number'].values.tolist()[0]
                true_ec_lst = uniecs.split('; ')
                for ec in true_ec_lst:
                    all_label.add(ec)
                true_label_dict[name] = true_ec_lst
                currpred_names.append(name)
                currpred_label.append(pred_label[i])
                currpred_probs.append(pred_probs[i])
                true_label.append(true_label_dict[name])
    #################### WITH ONLY UNIPROT ANNOTATION #################
    ###################################################################
    return true_label, all_label, currpred_names,currpred_label,currpred_probs
    
    
    # r_rm_index = rm_index.sort(reverse=True)
    # print('rm_index:',rm_index[:10],len(rm_index))
    # rm_index.sort(reverse=True)
    # # print('r_rm_index:',r_rm_index[:10])
    # for j in rm_index:
    #     del pred_names[j]
    #     del pred_label[j]
    #     del pred_probs[j]
    # # print('>>',len(pred_names),len(pred_label),len(pred_probs),len(all_label))
    # true_label = [true_label_dict[i] for i in pred_names]
    # return true_label, all_label, pred_names,pred_label,pred_probs

# def get_eval_metrics(pred_label, true_label, all_label):
#     mlb = MultiLabelBinarizer()
#     mlb.fit([list(all_label)])
#     n_test = len(pred_label)
#     pred_m = np.zeros((n_test, len(mlb.classes_)))
#     true_m = np.zeros((n_test, len(mlb.classes_)))
#     for i in range(n_test):
#         pred_m[i] = mlb.transform([pred_label[i]])
#         true_m[i] = mlb.transform([true_label[i]])
#     pre = precision_score(true_m, pred_m, average='weighted', zero_division=0)
#     rec = recall_score(true_m, pred_m, average='weighted')
#     f1 = f1_score(true_m, pred_m, average='weighted')
#     roc = roc_auc_score(true_m, pred_m, average='weighted')
#     acc = accuracy_score(true_m, pred_m)
#     return pre, rec, f1, roc, acc

def get_ec_pos_dict(mlb, true_label, pred_label):
    ec_list = []
    pos_list = []
    for i in range(len(true_label)):
        ec_list += list(mlb.inverse_transform(mlb.transform([true_label[i]]))[0])
        pos_list += list(np.nonzero(mlb.transform([true_label[i]]))[1])
    for i in range(len(pred_label)):
        ec_list += list(mlb.inverse_transform(mlb.transform([pred_label[i]]))[0])
        pos_list += list(np.nonzero(mlb.transform([pred_label[i]]))[1])
    label_pos_dict = {}
    for i in range(len(ec_list)):
        ec, pos = ec_list[i], pos_list[i]
        label_pos_dict[ec] = pos
        
    return label_pos_dict
def get_eval_metrics(pred_label, pred_probs, true_label, all_label):
    mlb = MultiLabelBinarizer()
    mlb.fit([list(all_label)])
    n_test = len(pred_label)
    pred_m = np.zeros((n_test, len(mlb.classes_)))
    true_m = np.zeros((n_test, len(mlb.classes_)))
    # for including probability
    pred_m_auc = np.zeros((n_test, len(mlb.classes_)))
    label_pos_dict = get_ec_pos_dict(mlb, true_label, pred_label)
    for i in range(n_test):
        pred_m[i] = mlb.transform([pred_label[i]])
        true_m[i] = mlb.transform([true_label[i]])
         # fill in probabilities for prediction
        labels, probs = pred_label[i], pred_probs[i]
        for label, prob in zip(labels, probs):
            if label in all_label:
                pos = label_pos_dict[label]
                pred_m_auc[i, pos] = prob
    pre = precision_score(true_m, pred_m, average='weighted', zero_division=0)
    rec = recall_score(true_m, pred_m, average='weighted')
    f1 = f1_score(true_m, pred_m, average='weighted')
    roc = roc_auc_score(true_m, pred_m_auc, average='weighted')
    acc = accuracy_score(true_m, pred_m)
    return pre, rec, f1, roc, acc

def pkl_getinfo(prediscoref,cutoff,all_label,refpred_names):
    dd = pd.read_pickle(prediscoref)
    # print(dd)
    pred_names, pred_probs, pred_label =[],[],[]
  
    pred_label_nocut=[]
    
    for name in refpred_names:
        pred_names.append(name)
        if name not in dd.keys():
            # probs=np.zeros(0)
            pred_probs.append([])
            pred_label.append([])
            pred_label_nocut.append([])
            continue
        val = dd[name]
        preds_ec_lst = []
        preds_ec_lst_nocut = []
        probs = np.zeros(len(dd.keys()))
        count = 0
        for ec,score in val.items():
            # ec_i = ec.split(":")[1].split("/")[0]
            ec_i = ec
            if ec_i not in all_label:
                # print('not in all ec ,ec_i:',ec_i)
                continue
            probs[count] = score
            count+=1
            preds_ec_lst_nocut.append(ec_i)
            if float(score) <= cutoff:
                preds_ec_lst.append(ec_i)     
   

        probs = (1 - np.exp(-1/probs)) / (1 + np.exp(-1/probs))
        probs = probs/np.sum(probs)

        pred_probs.append(probs)
        pred_label.append(preds_ec_lst)
        pred_label_nocut.append(preds_ec_lst_nocut)

    # for key, val in dd.items():
    #     key = key.split('=')[-1] ## remove after run tmp ecoli
    #     pred_names.append(key)
    #     preds_ec_lst = []
    #     preds_ec_lst_nocut = []
    #     probs = np.zeros(len(val.keys()))
    #     count = 0
    #     for ec,score in val.items():
    #         # ec_i = ec.split(":")[1].split("/")[0]
    #         ec_i = ec
    #         if ec_i not in all_label:
    #             continue
    #         probs[count] = score
    #         count+=1
    #         preds_ec_lst_nocut.append(ec_i)
    #         if float(score) <= cutoff:
    #             preds_ec_lst.append(ec_i)

    #     probs = (1 - np.exp(-1/probs)) / (1 + np.exp(-1/probs))
    #     probs = probs/np.sum(probs)

    #     pred_probs.append(probs)
    #     pred_label.append(preds_ec_lst)
    #     pred_label_nocut.append(preds_ec_lst_nocut)

    return pred_names,pred_label,pred_label_nocut,pred_probs



# 1	iLJ478	NC_000853.1	243274	Thermotoga maritima MSB8		243274	neg
# 2	iJN678	BA000022.2	1148	Synechocystis sp. PCC 6803		1148	neg
# 3	iJN1463	NC_002947.4	160488	Pseudomonas putida KT2440		160488	neg
# 4	iCN900	AM180355.1	272563	Clostridioides difficile 630		272563	pos
# 5	iHN637	NC_014328.1	748727	Clostridium ljungdahlii DSM 13528		748727	pos
# 6	iAF1260	NC_000913.3	83333	Escherichia coli str. K-12 substr. MG1655		511145 NOT WORKING BUT 83333 OK	neg
# 7	iAF987	CP000148.1	269799	Geobacter metallireducens GS-15		269799	neg
# 8	iCN718	NC_010410.1	509173	Acinetobacter baumannii AYE		509173	neg
# 9	iNJ661	NC_000962.3	83332	Mycobacterium tuberculosis H37Rv		83332	pos


# pref='/home/kexin/code/bigg/data/genomes/NC_000853.1'
# unif='/home/kexin/code/bigg/data/genomes/uniprot_243274.tsv'
# modelname = 'iLJ478'

# pref='/home/kexin/code/bigg/data/genomes/BA000022.2'
# unif='/home/kexin/code/bigg/data/genomes/uniprot_1148.tsv'
# modelname = 'iJN678'

# pref='/home/kexin/code/bigg/data/genomes/NC_002947.4'
# unif='/home/kexin/code/bigg/data/genomes/uniprot_160488.tsv'
# modelname = 'iJN1463'

# pref='/home/kexin/code/bigg/data/genomes/AM180355.1'
# unif='/home/kexin/code/bigg/data/genomes/uniprot_272563.tsv'
# modelname = 'iCN900'

# pref='/home/kexin/code/bigg/data/genomes/NC_014328.1'
# unif='/home/kexin/code/bigg/data/genomes/uniprot_748727.tsv'
pref='/ibex/user/niuk0a/CLEAN/app/results/inputs/NC_014328.1'
unif='/ibex/user/niuk0a/funcarve/cobra/uniprot_748727.tsv'
modelname = 'iHN637'

# pref='/home/kexin/code/bigg/data/genomes/NC_000913.3'
# unif='/home/kexin/code/bigg/data/genomes/uniprot_83333.tsv'
# modelname = 'iAF1260'

# pref='/home/kexin/code/bigg/data/genomes/CP000148.1'
# unif='/home/kexin/code/bigg/data/genomes/uniprot_269799.tsv'
# modelname = 'iAF987'

# pref='/home/kexin/code/bigg/data/genomes/NC_010410.1'
# unif='/home/kexin/code/bigg/data/genomes/uniprot_509173.tsv'
# modelname = 'iCN718'

# pref='/home/kexin/code/bigg/data/genomes/NC_000962.3'
# unif='/home/kexin/code/bigg/data/genomes/uniprot_83332.tsv'
# modelname = 'iNJ661'

print('genome_name:',pref.split('/')[-1])
result =[]

pred_label, pred_names, pred_all_label = get_pred_labels(pref)
pred_probs = get_pred_probs(pref)
# print('ori',len(pred_label),len(pred_names),len(pred_all_label),len(pred_probs))
true_label, all_label, pred_names, pred_label, pred_probs = get_true_labels(unif,pref,pred_names,pred_label,pred_probs)
# print('ori',len(true_label),len(all_label),len(pred_names),len(pred_label),len(pred_probs))
# 
pre,rec,f1,roc,acc = get_eval_metrics(pred_label, pred_probs, true_label, all_label)
tmpresult = 'CLEAN|'+str(pre)+'|'+str(rec)+'|'+str(f1)+'|'+str(roc)+'|'+str(acc)+'|'+str(len(pred_names))+'|'+str(len(all_label))
# print('ori:',pre,rec,f1,roc,acc,'|',len(pred_names),len(all_label))
result.append(tmpresult)
# for i in range(1,10):
for i in range(1,5):
    # 
    # /ibex/user/niuk0a/funcarve/cobra/iHN637_R_newpredscore_7.pkl
    newprediscoref = '/ibex/user/niuk0a/funcarve/cobra/'+str(modelname) +'_R_newpredscore_'+str(i)+'.pkl'
    # newprediscoref = '/home/kexin/code/bigg/data/genomes/results/' + str(modelname) +'_R_newpredscore_'+str(i)+'.pkl'
    interpred_names,interpred_label_cutoff,interpred_label,interpred_probs = pkl_getinfo(newprediscoref,0.2,all_label,pred_names)
    # print('inter',len(true_label),len(all_label),len(interpred_names),len(interpred_label),len(interpred_probs))
    inter_pre,inter_rec,inter_f1,inter_roc,inter_acc = get_eval_metrics(interpred_label, interpred_probs, true_label, all_label) ## INTER
    tmpresult = 'Interaction '+str(i)+'|'+str(inter_pre)+'|'+str(inter_rec)+'|'+str(inter_f1)+'|'+str(inter_roc)+'|'+str(inter_acc)+'|'+str(len(interpred_names))+'|'+str(len(all_label))
    # print('inter'+str(i)+':',inter_pre,inter_rec,inter_f1,inter_roc,inter_acc,'|',len(pred_names),len(all_label))
    result.append(tmpresult)

print('genome_name|',pref.split('/')[-1],'|',modelname,sep='')
for r in result:
    print(r)
    



