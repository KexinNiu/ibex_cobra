from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.metrics import precision_score, recall_score, \
    roc_auc_score, accuracy_score, f1_score, average_precision_score
import pandas as pd
import csv
import numpy as np
import os
from Bio import SeqIO
import torch

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
    print(f'>pred_label:{len(pred_label)},pred_names:{len(pred_names)},pred_all_label:{len(pred_all_label)}')
    return pred_label, pred_names,pred_all_label

def get_pred_probs(out_filename, pred_type="_maxsep"):
    file_name = out_filename+pred_type
    result = open(file_name+'.csv', 'r')
    csvreader = csv.reader(result, delimiter=',')
    pred_probs = []
    # for row in csvreader:
    #     preds_ec_lst = []
    #     preds_with_dist = row[1:]
    #     probs = torch.zeros(len(preds_with_dist))
    #     # probs = np.zeros(len(preds_with_dist))
    #     count = 0
    #     for pred_ec_dist in preds_with_dist:
    #         # get EC number 3.5.2.6 from EC:3.5.2.6/10.8359
    #         ec_i = float(pred_ec_dist.split(":")[1].split("/")[1])
    #         probs[count] = ec_i
    #         #preds_ec_lst.append(probs)
    #         count += 1
    #     # sigmoid of the negative distances 
    #     # probs = (1 - torch.exp(-1/probs)) / (1 + torch.exp(-1/probs))
    #     # print('probs:',probs,end='->',flush=True)
    #     probs = (1 - torch.exp(-1/probs)) / (1 + torch.exp(-1/probs))
    #     probs = probs/torch.sum(probs)
    #     # print('probs:',probs,flush=True)

    #     pred_probs.append(probs)
    for row in csvreader:
        preds_with_dist = row[1:]
        probs = []
        for pred_ec_dist in preds_with_dist:
            # 提取距离信息
            dist = float(pred_ec_dist.split(":")[1].split("/")[1])
            probs.append(-dist)  # 取负数作为logits，用于softmax
        
        # 使用softmax归一化以生成概率
        probs = torch.tensor(probs)
        probs = torch.softmax(probs, dim=0).numpy()
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

# def get_pred_labels_prc(out_filename, cutoff, pred_type="_maxsep"):
#     file_name = out_filename+pred_type
#     result = open(file_name+'.csv', 'r')
#     csvreader = csv.reader(result, delimiter=',')
#     pred_label = []
#     pred_all_label=set()
#     for row in csvreader:
#         preds_ec_lst = []
#         preds_with_dist = row[1:]
#         pred_names.append(row[0])

#         for pred_ec_dist in preds_with_dist:
#             # get EC number 3.5.2.6 from EC:3.5.2.6/10.8359
#             ec_i = pred_ec_dist.split(":")[1].split("/")[0]
#             if int(pred_ec_dist.split(":")[1].split("/")[1]) <= cutoff:
#                 preds_ec_lst.append(ec_i)
#                 pred_all_label.add(ec_i)
#         pred_label.append(preds_ec_lst)
#     return pred_label ,pred_names,pred_all_label

def get_true_labels(file_name,pref,pred_names,pred_label,pred_probs,pred_all_label):
    
    # gbkf = pref+'.gb'
    gbkf = pref[:-3] + '.gb'
    print('gbkf:',gbkf)
    if not os.path.exists(gbkf):
        print('no gbkf:!',gbkf)
    r_df,df,lt2ec,lt2oldlt, oldlt2lt = read_unif_genebk(file_name,gbkf)
    all_label = set()
    true_label_dict = {}
    rm_index = []
    goinfo={}
    true_label = []
    
    currpred_names = []
    currpred_label = []
    currpred_probs = []
    # print('>>',len(pred_names),len(pred_label),len(pred_probs),len(all_label))
    print(f'get_true_labels:len(pred_label):{len(pred_label)},len(pred_names):{len(pred_names)},len(pred_all_label):{len(pred_all_label)},len(pred_probs):{len(pred_probs)}')
    if not lt2ec == {}:
    #########################################################################
    #################### WITH GENEBANK + UNIPROT ANNOTATION #################
        print('lt2ec!=empty')
        # print('lt2oldlt:',lt2oldlt)
        for i in range(0,len(pred_names)):
            name = pred_names[i]
            try:
                oldname = lt2oldlt[name]
                # print('oldname:',oldname)
            except:
                # print('no oldname:',name)
                #stop here
                pass
                # print('skip',name)
                
            remove = True
            for on in oldname:
                ## for each gene names, split by ' ' if there is on , then return unip
                unip  = r_df[r_df['Gene Names'].str.contains(on+' ', na=False) | r_df['Gene Names'].str.endswith(on, na=False)]
                # unip = r_df[r_df['Gene Names'].str.contains(on, na=False)]
                if on in lt2ec.keys():
                    if len(lt2ec[name]) > 0:
                        remove = False
                        true_ec_lst = lt2ec[name]
                        for ec in true_ec_lst:
                            all_label.add(ec)
                        true_label_dict[name] = true_ec_lst
                elif len(unip) != 0:
                    if unip['EC number'].isnull().values.any():
                        # print('unip ec number null')
                        continue
                    else:
                        remove = False
                        uniecs = unip['EC number'].values.tolist()[0]
                        true_ec_lst = uniecs.split('; ')
                        for ec in true_ec_lst:
                            all_label.add(ec)
                        true_label_dict[name] = true_ec_lst
                # if len(unip) == 0:
                #     print('unip len 0') 
                #     continue
                # elif unip['EC number'].isnull().values.any():
                #     print('unip ec number null')
                #     continue
                # elif not unip['EC number'].isnull().values.any():
                #     remove = False
                #     uniecs = unip['EC number'].values.tolist()[0]
                #     true_ec_lst = uniecs.split('; ')
                #     for ec in true_ec_lst:
                #         all_label.add(ec)
                #     true_label_dict[name] = true_ec_lst
                # elif len(lt2ec[on]) > 0:
                #     remove = False
                #     true_ec_lst = lt2ec[on]
                #     for ec in true_ec_lst:
                #         all_label.add(ec)
                #     true_label_dict[name] = true_ec_lst
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
        print(' WITH ONLY UNIPROT ANNOTATION ')        
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
        smallest_10_dist_df = pd.Series(val).nsmallest(20)
        dist_lst = list(smallest_10_dist_df)
        max_sep_i = maximum_separation(dist_lst, True, True)
        max_sep_ec = [smallest_10_dist_df.index[i] for i in range(max_sep_i+1)]

        for ec,score in val.items():
            # ec_i = ec.split(":")[1].split("/")[0]
            ec_i = ec
            if ec_i not in all_label:
                # print('not in all ec ,ec_i:',ec_i)
                continue
            if ec_i not in max_sep_ec:
                continue
            probs[count] = score
            count+=1
            preds_ec_lst_nocut.append(ec_i)
            if float(score) <= cutoff:
                preds_ec_lst.append(ec_i)     
   
        epsilon = 1e-10
        probs = (1 - np.exp(-1/(probs + epsilon))) / (1 + np.exp(-1/(probs + epsilon)))
        # probs = (1 - np.exp(-1/probs)) / (1 + np.exp(-1/probs))
        probs = probs/np.sum(probs)

        pred_probs.append(probs)
        pred_label.append(preds_ec_lst)
        pred_label_nocut.append(preds_ec_lst_nocut)

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
## allec reward =0.1
# 		genome 	ORGID	ORGNASIM	
# t1	iLJ478	NC_000853.1	243274	Thermotoga maritima MSB8	neg
# t2	iJN1463	NC_002947.4	160488	Pseudomonas putida KT2440	neg
# t2	iCN900	AM180355.1	272563	Clostridioides difficile 630	pos
# t4	iAF1260	NC_000913.3	83333	Escherichia coli str. K-12 substr. MG1655	neg
# t4	iAF987	CP000148.1	269799	Geobacter metallireducens GS-15	neg
# t5	iNJ661	NC_000962.3	83332	Mycobacterium tuberculosis H37Rv	pos
# t4	iYO844	AL009126.3	224308	Bacillus subtilis subsp. subtilis str. 168	pos

# pref='/ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1'
# unif='/ibex/user/niuk0a/funcarve/cobra//uniprot_243274.tsv'
# modelname = 'iLJ478'

# pref='/ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2'
# unif='/ibex/user/niuk0a/funcarve/cobra//uniprot_160488.tsv'
# modelname = 'iJN1463'


# pref='/ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1'
# unif='/ibex/user/niuk0a/funcarve/cobra//uniprot_272563.tsv'
# modelname = 'iCN900'

# pref='/ibex/user/niuk0a/CLEAN/app/results/inputs/NC_014328.1'
# unif='/ibex/user/niuk0a/funcarve/cobra//uniprot_748727.tsv'
# pref='/ibex/user/niuk0a/CLEAN/app/results/inputs/NC_014328.1'
# unif='/ibex/user/niuk0a/funcarve/cobra/uniprot_748727.tsv'
# modelname = 'iHN637'

# pref='/ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3'
# unif='/ibex/user/niuk0a/funcarve/cobra//uniprot_83333.tsv'
# modelname = 'iAF1260'

# pref='/ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1'
# unif='/ibex/user/niuk0a/funcarve/cobra//uniprot_269799.tsv'
# modelname = 'iAF987'

# pref='/ibex/user/niuk0a/CLEAN/app/results/inputs/NC_010410.1'
# unif='/ibex/user/niuk0a/funcarve/cobra//uniprot_509173.tsv'
# modelname = 'iCN718'

# pref='/ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3'
# unif='/ibex/user/niuk0a/funcarve/cobra//uniprot_83332.tsv'
# modelname = 'iNJ661'

def eval(pref,unif,modelname,thrname,evaltype,of,resultmatrix,all_label,pred_names,df):

    print('genome_name:',pref.split('/')[-1])
    if 'p' in thrname:
        thrval = float('0.'+thrname.split('p')[-1][-1])
    else:
        thrval = float(thrname)
    

    for i in range(1,10):
        if len(df) >0:
            if df[(df['modelname']==modelname) & (df['type']==evaltype) & (df['thr']==thrval) & (df['iter_round']==i)].shape[0] > 0:
                print('already done skip:',modelname,evaltype,thr)
                continue
        newprediscoref = f'/ibex/user/niuk0a/funcarve/cobra/{modelname}_{evaltype}_{thrname}_t6_newpredscore_{i}.pkl'
        if not os.path.exists(newprediscoref):
            continue
        interpred_names,interpred_label_cutoff,interpred_label,interpred_probs = pkl_getinfo(newprediscoref,0.2,all_label,pred_names)
        # print('inter',len(true_label),len(all_label),len(interpred_names),len(interpred_label),len(interpred_probs))
        try:
            inter_pre,inter_rec,inter_f1,inter_roc,inter_acc = get_eval_metrics(interpred_label, interpred_probs, true_label, all_label)
        # tmpresult = 'Interaction '+str(i)+'|'+str(inter_pre)+'|'+str(inter_rec)+'|'+str(inter_f1)+'|'+str(inter_roc)+'|'+str(inter_acc)+'|'+str(len(interpred_names))+'|'+str(len(all_label))
        # print('inter'+str(i)+':',inter_pre,inter_rec,inter_f1,inter_roc,inter_acc,'|',len(pred_names),len(all_label))
            
            resultmatrix['modelname'].append(modelname)
            resultmatrix['type'].append(evaltype)
            resultmatrix['iter_round'].append(i)
            resultmatrix['thr'].append(thrval)
            resultmatrix['pre'].append(inter_pre)
            resultmatrix['rec'].append(inter_rec)
            resultmatrix['f1'].append(inter_f1)
            resultmatrix['roc'].append(inter_roc)
            resultmatrix['acc'].append(inter_acc)
            resultmatrix['pred_len'].append(len(interpred_names))
            resultmatrix['all_len'].append(len(all_label))
        except:
            continue
    return resultmatrix
    # print('genome_name|',pref.split('/')[-1],'|',modelname,'|',evaltype,sep='',file=of)

    # for r in result:
    #     print(r,file=of,flush=True)
# print('genome_name:',pref.split('/')[-1])
# result =[]

# pred_label, pred_names, pred_all_label = get_pred_labels(pref)
# pred_probs = get_pred_probs(pref)
# # print('ori',len(pred_label),len(pred_names),len(pred_all_label),len(pred_probs))
# true_label, all_label, pred_names, pred_label, pred_probs = get_true_labels(unif,pref,pred_names,pred_label,pred_probs)
# # print('ori',len(true_label),len(all_label),len(pred_names),len(pred_label),len(pred_probs))
# # 
# pre,rec,f1,roc,acc = get_eval_metrics(pred_label, pred_probs, true_label, all_label)
# tmpresult = 'CLEAN|'+str(pre)+'|'+str(rec)+'|'+str(f1)+'|'+str(roc)+'|'+str(acc)+'|'+str(len(pred_names))+'|'+str(len(all_label))
# # print('ori:',pre,rec,f1,roc,acc,'|',len(pred_names),len(all_label))
# result.append(tmpresult)
# # for i in range(1,10):
# for i in range(1,10):
#     # /ibex/user/niuk0a/funcarve/cobra/iYO844_allec_t8_newpredscore_10.pkl
#     # /ibex/user/niuk0a/funcarve/cobra/iHN637_R_newpredscore_7.pkl
#     newprediscoref = '/ibex/user/niuk0a/funcarve/cobra/'+str(modelname) +'_R_newpredscore_'+str(i)+'.pkl'
#     # newprediscoref = '/home/kexin/code/bigg/data/genomes/results/' + str(modelname) +'_R_newpredscore_'+str(i)+'.pkl'
#     interpred_names,interpred_label_cutoff,interpred_label,interpred_probs = pkl_getinfo(newprediscoref,0.2,all_label,pred_names)
#     # print('inter',len(true_label),len(all_label),len(interpred_names),len(interpred_label),len(interpred_probs))
#     inter_pre,inter_rec,inter_f1,inter_roc,inter_acc = get_eval_metrics(interpred_label, interpred_probs, true_label, all_label) ## INTER
#     tmpresult = 'Interaction '+str(i)+'|'+str(inter_pre)+'|'+str(inter_rec)+'|'+str(inter_f1)+'|'+str(inter_roc)+'|'+str(inter_acc)+'|'+str(len(interpred_names))+'|'+str(len(all_label))
#     # print('inter'+str(i)+':',inter_pre,inter_rec,inter_f1,inter_roc,inter_acc,'|',len(pred_names),len(all_label))
#     result.append(tmpresult)

# print('genome_name|',pref.split('/')[-1],'|',modelname,sep='')
# for r in result:
#     print(r)
    



# 		genome 	ORGID	ORGNASIM	
# t1	iLJ478	NC_000853.1	243274	Thermotoga maritima MSB8	neg
# t2	iJN1463	NC_002947.4	160488	Pseudomonas putida KT2440	neg
# t2	iCN900	AM180355.1	272563	Clostridioides difficile 630	pos
# t4	iAF1260	NC_000913.3	83333	Escherichia coli str. K-12 substr. MG1655	neg
# t4	iAF987	CP000148.1	269799	Geobacter metallireducens GS-15	neg
# t5	iNJ661	NC_000962.3	83332	Mycobacterium tuberculosis H37Rv	pos
# t4	iYO844	AL009126.3	224308	Bacillus subtilis subsp. subtilis str. 168	pos

# pref='/ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1'
# unif='/ibex/user/niuk0a/funcarve/cobra/uniprot_243274.tsv'
# modelname = 'iLJ478'

# pref='/ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2'
# unif='/ibex/user/niuk0a/funcarve/cobra/uniprot_160488.tsv'
# modelname = 'iJN1463'

prefs =[
    '/ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1',
    '/ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2',
    '/ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2',
    '/ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3_t4',
    '/ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_t4',
    '/ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_t5',
    '/ibex/user/niuk0a/CLEAN/app/results/inputs/AL009126.3_t4'
]
unifs = [   
    '/ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_243274.tsv',
    '/ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_160488.tsv',
    '/ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_272563.tsv',
    '/ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_83333.tsv',
    '/ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_269799.tsv',
    '/ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_83332.tsv',
    '/ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_224308.tsv'
]

modelnames = [
    'iLJ478',
    'iJN1463',
    'iCN900',
    'iAF1260',
    'iAF987',
    'iNJ661',
    'iYO844'
]
# outf='/ibex/user/niuk0a/funcarve/cobra/eval_clean_all_7.out'
resultmatrix = {
    'modelname':[],
    'type':[],
    'iter_round':[],
    'thr':[],
    'pre':[],
    'rec':[],
    'f1':[],
    'roc':[],
    'acc':[],
    'pred_len':[],
    'all_len':[]
}
# import argparse
# parser = argparse.ArgumentParser(description='Generate genome-scale metabolic network reconstruction from KEGG BLAST hits.')
# parser.add_argument('--evaltype', default='allec', type=str, help='allec,flux,fluxblock,block')
# evaltype = parser.parse_args().evaltype
evallist=['allec','flux','fluxblock','block']
# evaltype='flux'
# allec flux fluxblock block
# outf = f'/ibex/user/niuk0a/funcarve/cobra/eval_clean_all_7_{evaltype}.out'
# of = open(outf,'w')
file_path = '/ibex/user/niuk0a/funcarve/cobra/eval_clean_7v3.csv'
file_path = '/ibex/user/niuk0a/funcarve/cobra/data/result/eval_clean_7v4.csv'
if not os.path.exists(file_path):
    df = pd.DataFrame()
else:
    df = pd.read_csv(file_path)
    

print('df:',df.head())

for i in range(0,len(prefs)):
    pref = prefs[i]
    unif = unifs[i]
    modelname = modelnames[i]

    pred_label, pred_names, pred_all_label = get_pred_labels(pref)
    pred_probs = get_pred_probs(pref)
    true_label, all_label, pred_names, pred_label, pred_probs = get_true_labels(unif,pref,pred_names,pred_label,pred_probs,pred_all_label)
    print('ori',len(true_label),len(all_label),len(pred_names),len(pred_label),len(pred_probs))
    ## check if we already have the result in df
    if len(df) > 0:
        if df[(df['modelname']==modelname) & (df['type']=='CLEAN') & (df['thr']==0) & (df['iter_round']==0)].shape[0] > 0:
            print('already done skip:',modelname)
    else:
        pre,rec,f1,roc,acc = get_eval_metrics(pred_label, pred_probs, true_label, all_label)
    # tmpresult = 'CLEAN|'+str(pre)+'|'+str(rec)+'|'+str(f1)+'|'+str(roc)+'|'+str(acc)+'|'+str(len(pred_names))+'|'+str(len(all_label))
        resultmatrix['modelname'].append(modelname)
        resultmatrix['type'].append('CLEAN')
        resultmatrix['iter_round'].append(0)
        resultmatrix['thr'].append(0)
        resultmatrix['pre'].append(pre)
        resultmatrix['rec'].append(rec)
        resultmatrix['f1'].append(f1)
        resultmatrix['roc'].append(roc)
        resultmatrix['acc'].append(acc)
        resultmatrix['pred_len'].append(len(pred_names))
        resultmatrix['all_len'].append(len(all_label))

    # thrlist = [0.1,0.3,0.5,0.7,1.0,0.05, 0.2,2.0, 3.0,5.0]
    # thrlist = [0.1,0.3,0.5,0.7,1.0]rewards=(0.01 1.5 2.0 3.0 5.0 )
    # thrlist = [0.1,0.3,0.5,0.7,1.0,0.01, 1.5, 2.0 ,3.0 ,5.0]
    # thrlist = [0.01, 2.0 ,3.0 ,5.0]
    thrlist = [ 0.1,0.4,0.8,1.0 ,3.0]
    
    for thr in thrlist:
        if thr < 1.0:
            thrname = 'p'+str(thr).replace('.','')
        else:
            thrname = str(thr).split('.')[0]
            if str(thr).split('.')[1] != '0':
                thrname = thrname+'p'+str(thr).split('.')[1]
        print('thrname:',thrname)
        for evaltype in evallist:
            ## check if we already have the result in df
            if len(df) > 0:
                if df[(df['modelname']==modelname) & (df['type']==evaltype) & (df['thr']==thr) & (df['iter_round']==0)].shape[0] > 0:
                    print('already done skip:',modelname,evaltype,thr)
                
                    continue
            else:
                print('start:',modelname,evaltype,thr)
            # resultmatrix = eval(pref,unif,modelname,thrname,evaltype,of,resultmatrix,all_label,pred_names)
                try:
                    resultmatrix = eval(pref,unif,modelname,thrname,evaltype,'of',resultmatrix,all_label,pred_names,df)
                except:
                    print('error:',modelname,evaltype,thr)
                    continue
                print('*'*10)
   

new_df = pd.DataFrame(resultmatrix)
if len(df) == 0:
    new_df.to_csv(file_path, index=False)
else:
    df = pd.concat([df, new_df], ignore_index=True)
    df.to_csv(file_path, index=False)


# new_df.to_csv(file_path, index=False)
# df.to_csv(f'/ibex/user/niuk0a/funcarve/cobra/eval_clean_7v3.csv',index=False)