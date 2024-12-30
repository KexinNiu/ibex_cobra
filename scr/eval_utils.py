import pickle
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

def get_true_labels(file_name,pref,pred_names,pred_label,pred_probs,pred_all_label):
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
    print(f'get_true_labels:len(pred_label):{len(pred_label)},len(pred_names):{len(pred_names)},len(pred_all_label):{len(pred_all_label)},len(pred_probs):{len(pred_probs)}')
    print('pred names:',pred_names[:5])
    if not lt2ec == {}:
    #########################################################################
    #################### WITH GENEBANK + UNIPROT ANNOTATION #################
        print('WITH GENEBANK + UNIPROT ANNOTATION')
        for i in range(0,len(pred_names)):
            name = pred_names[i]
            try:
                oldname = lt2oldlt[name]
                # print('oldname:',oldname)

            except:
                pass
                
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
                            if '-' in ec:
                                continue
                            all_label.add(ec)
                        true_ec_lst = [ec for ec in true_ec_lst if '-' not in ec]
                        if len(true_ec_lst) > 0:
                            true_label_dict[name] = true_ec_lst
                            # print('name with label found:',name,true_ec_lst)
                        else:
                            remove = True
                            
              
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
       
def read_unif(f):
    df = pd.read_csv(f,sep='\t')
    # only take where Reviewed == reviewed
    r_df = df[df['Reviewed']=='reviewed'] 
    #  take columns ['Entry','Protein names','Gene names','Organism','EC number','Status']
    r_df = r_df[['Entry','Gene Names','EC number','Gene Ontology (GO)']]
    return r_df,df

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

# def update_differ(oridf,differdf):
#     dd = oridf.copy()
#     for key in differdf.keys():
#         dd[key] = differdf[key]
#     return dd

def update_differ(oridf,differdf):
    dd = oridf.copy()
    for pr in differdf.keys():
        for ec in differdf[pr].keys():
            if pr not in dd.keys():
                dd[pr] = {ec:differdf[pr][ec]}
            else:
                dd[pr][ec] = differdf[pr][ec]
    return dd

def pr2ec_to_ec2pr(allpr2ec): ## predscore
    allec2pr = {}
    for pr in allpr2ec.keys():
        for ec in allpr2ec[pr].keys():
            try:
                allec2pr[ec].update({pr:allpr2ec[pr][ec]})
            except:
                allec2pr[ec] = {pr:allpr2ec[pr][ec]}
    return allec2pr
            
def maximum_separation(dist_lst, first_grad, use_max_grad):
    #cite from CLEAN
    opt = 0 if first_grad else -1
    gamma = np.append(dist_lst[1:], np.repeat(dist_lst[-1], 10))
    sep_lst = np.abs(dist_lst - np.mean(gamma))
    sep_grad = np.abs(sep_lst[:-1]-sep_lst[1:])
    if use_max_grad:
        max_sep_i = np.argmax(sep_grad)
    else:
        large_grads = np.where(sep_grad > np.mean(sep_grad))
        max_sep_i = large_grads[-1][opt]
    if max_sep_i >= 5:
        max_sep_i = 0
    return max_sep_i

def infer_confidence_gmm(distance, gmm_lst):
    confidence = []
    for j in range(len(gmm_lst)):
        main_GMM = gmm_lst[j]
        a, b = main_GMM.means_
        true_model_index = 0 if a[0] < b[0] else 1
        certainty = main_GMM.predict_proba([[distance]])[0][true_model_index]
        confidence.append(certainty)
    return np.mean(confidence)

def write_max_sep_choices(df, csv_name, first_grad=True, use_max_grad=False, gmm = None):
    out_file = open(csv_name + '_maxsep.csv', 'w', newline='')
    csvwriter = csv.writer(out_file, delimiter=',')
    all_test_EC = set()
    for col in df.columns:
        ec = []
        smallest_10_dist_df = df[col].nsmallest(10)
        dist_lst = list(smallest_10_dist_df)
        max_sep_i = maximum_separation(dist_lst, first_grad, use_max_grad)
        for i in range(max_sep_i+1):
            EC_i = smallest_10_dist_df.index[i]
            dist_i = smallest_10_dist_df[i]
            if gmm != None:
                gmm_lst = pickle.load(open(gmm, 'rb'))
                dist_i = infer_confidence_gmm(dist_i, gmm_lst)
            dist_str = "{:.4f}".format(dist_i)
            all_test_EC.add(EC_i)
            ec.append('EC:' + str(EC_i) + '/' + dist_str)
        ec.insert(0, col)
        csvwriter.writerow(ec)
    return

def eval_df(dd,all_label,refpred_names,cutoff=5):
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
            ec_i = ec
            if ec_i not in all_label:
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

def get_ancester(eclist):
    allecs = set()
    for ec in eclist:
        digs =  ec.split('.')
        for i in range(4):
            allecs.add('.'.join(digs[:i+1]))
    return allecs

def get_metric(ec_pred,ec_true):
    ec_pred = get_ancester(ec_pred)
    ec_true = get_ancester(ec_true)
    union_len = len(ec_pred.union(ec_true))
    if union_len == 0:
        print(ec_pred,ec_true)
    true_len = len(ec_true)
    pred_len = len(ec_pred)
    acc = len(ec_pred.intersection(ec_true)) / union_len if union_len != 0 else 0
    recall = len(ec_pred.intersection(ec_true)) / true_len if true_len != 0 else 0
    precision = len(ec_pred.intersection(ec_true)) / pred_len if pred_len != 0 else 0
    f1_score = 2 * recall * precision / (recall + precision) if (recall + precision) != 0 else 0
    return acc, recall, precision, f1_score


