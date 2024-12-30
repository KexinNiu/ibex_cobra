import pickle
import pandas as pd
import numpy as np
import os
import csv
import torch
from sklearn.metrics import (
    precision_score, recall_score, f1_score, roc_auc_score, average_precision_score, 
    mean_squared_error, log_loss, brier_score_loss
)
folder = '/ibex/user/niuk0a/funcarve/cobra/data/tmp'


f='/ibex/user/niuk0a/funcarve/cobra/split70_bactaxo_test2.csv'
df = pd.read_csv(f,sep='\t')
prs = df['Entry'].values.tolist()  # 12408
print('prs',prs[:10],len(prs))  # 12408

def get_all_data(prs, folder):
    data = {}
    prcount = {pr:0 for pr in prs}
    for file in os.listdir(folder):
        if file.startswith('UP') and file.endswith('.pkl') and 'ORIpredscore' in file:
            index = file.split('_')[0]
            print('index',index)
            with open(os.path.join(folder, file), 'rb') as f:
                currentdata = pickle.load(f)
                # format keys split by '|'keeps 1
                for pr in list(currentdata.keys()):
                    newpr = pr.split('|')[1]
                    if newpr in prs:
                        currentdata[newpr] = currentdata.pop(pr)
                    else:
                        currentdata.pop(pr)                
            
            for i in range(1,4):
                difff = file.replace('ORIpredscore_0','differ_predscore_'+str(i))
                with open(os.path.join(folder, difff), 'rb') as f:
                    dd = pickle.load(f)
                    for pr,info in dd.items():
                        pr = pr.split('|')[1]
                        if pr in prs:
                            for ec,val in info.items():
                                currentdata[pr][ec]=val
            for key in currentdata.keys():
                prcount[key] += 1
            data[index] = currentdata
    return data,prcount
def get_true_labels_mapname(pred_names,file_name):
    result = open(file_name+'.csv', 'r')
    csvreader = csv.reader(result, delimiter='\t')
    all_label = set()
    true_label_dict = {}
    header = True
    count = 0
    for row in csvreader:
        # don't read the header
        if header is False:
            if row[0] not in pred_names:
                continue
            count += 1
            true_ec_lst = row[1].split(';')
            true_label_dict[row[0]] = true_ec_lst
            for ec in true_ec_lst:
                all_label.add(ec)
        if header:
            header = False
    true_label = [true_label_dict[i] for i in pred_names]
    return true_label, all_label

t2_cleandf = pd.read_pickle('/ibex/user/niuk0a/CLEAN/app/results/inputs/split70_bactaxo_test2_maxsep_df.pkl')
t2_iterdf = pd.read_pickle('/ibex/user/niuk0a/funcarve/cobra/data/tmp/up_combine_data.pkl')
t2_iterdf = t2_iterdf.T


cleancols = t2_cleandf.columns.to_list()
itercols = t2_iterdf.columns.to_list()
print(t2_cleandf.shape)
print(t2_cleandf.head())
exit()
multicols = [i for i in itercols if i in prs]
print('multicols:',len(multicols))

t2_cleandf = t2_cleandf[multicols]
t2_iterdf = t2_iterdf[multicols]

pred_names = t2_iterdf.columns.to_list()
print('pred_names:',len(pred_names),pred_names[:10])
#### scale the data

# mask = (t2_iterdf != t2_cleandf)
# print('type:',type(mask),mask.shape)
# # print('mask:',mask)

# t2_iterdf[mask] = t2_iterdf[mask] -0.5
# print ('t2_iterdf:',t2_iterdf.shape,'head:\n',t2_iterdf.head(),sep='')


########################



t2_cleandf = t2_cleandf[sorted(t2_cleandf.columns)]
t2_iterdf = t2_iterdf[sorted(t2_iterdf.columns)]
pred_names = t2_iterdf.columns.to_list()


from eval_utils import get_eval_metrics, get_pred_labels, get_pred_probs, write_max_sep_choices, infer_confidence_gmm

def get_pred_labels(out_filename, pred_type="_maxsep"):
    file_name = out_filename+pred_type
    result = open(file_name+'.csv', 'r')
    csvreader = csv.reader(result, delimiter=',')
    pred_label = []
    pred_names = []
    for row in csvreader:
        preds_ec_lst = []
        preds_with_dist = row[1:]
        pred_names.append(row[0])
        for pred_ec_dist in preds_with_dist:
            # get EC number 3.5.2.6 from EC:3.5.2.6/10.8359
            ec_i = pred_ec_dist.split(":")[1].split("/")[0]
            preds_ec_lst.append(ec_i)
        pred_label.append(preds_ec_lst)
    return pred_names,pred_label

true_label,all_label = get_true_labels_mapname(pred_names,'/ibex/user/niuk0a/funcarve/cobra/split70_bactaxo_test2')

print('true_label:',len(true_label),'all_label',len(all_label))
print('true_label:',true_label[:10])    
all_label = list(all_label)
cleanacc,cleanrecall,cleanpre,cleanf1 = 0,0,0,0
iteracc,iterrecall,iterpre,iterf1 = 0,0,0,0
from eval_utils import get_metric,get_ancester
out_filename = '/ibex/user/niuk0a/funcarve/cobra/data/result/split70t2_ITERt2v4'
it_names,pred_label = get_pred_labels(out_filename, pred_type='_maxsep')
print('it_names:',len(it_names),it_names[:10])
for i, pred in enumerate(pred_names):
    acc,recall,pre,f1 = get_metric(pred_label[it_names.index(pred)],true_label[i])
    iteracc += acc
    iterrecall += recall
    iterpre += pre
    iterf1 += f1

  
out_filename ='/ibex/user/niuk0a/funcarve/cobra/data/result/split70t2_CLEANt2_subset'
c_names,pred_label = get_pred_labels(out_filename, pred_type='_maxsep')
print('pred_names:',len(pred_names),pred_names[:10])
for i, pred in enumerate(pred_names):
    acc,recall,pre,f1 = get_metric(pred_label[c_names.index(pred)],true_label[i])
    cleanacc += acc
    cleanrecall += recall
    cleanpre += pre
    cleanf1 += f1
    
print('clean\t',cleanacc/len(pred_names),cleanrecall/len(pred_names),cleanpre/len(pred_names),cleanf1/len(pred_names))


print('type\tacc\trecall\tpre\tf1')
print('clean\t',cleanacc/len(pred_names),cleanrecall/len(pred_names),cleanpre/len(pred_names),cleanf1/len(pred_names))
print('iter\t',iteracc/len(pred_names),iterrecall/len(pred_names),iterpre/len(pred_names),iterf1/len(pred_names))
# type    acc     recall  pre     f1
# clean    0.9819451273052388 0.987147052572953 0.9884412640299071 0.9867052997817197
# iter     0.9816604243197049 0.9869739223790472 0.9883066072124248 0.9865244749125291
exit()

# label = torch.zeros(len(pred_names), len(all_label))

# for i, pred in enumerate(pred_names):
#     # print('pred protein:',pred)
#     for ec in true_label[i]:
#         j = all_label.index(ec)
#         label[i, j] = 1

# print('label:',label.shape)
# print('label:',label[0][:10])



# gmm = '/ibex/user/niuk0a/CLEAN/app/data/pretrained/gmm_t2.pkl'
# gmm_lst = pickle.load(open(gmm, 'rb'))


########## clean data # t2_cleandf
# print('t2_cleandf:',t2_cleandf.shape,'head:\n',t2_cleandf.head(),sep='')
# prelabel = torch.zeros(len(pred_names), len(all_label))

# for i, pred in enumerate(pred_names):
#     subdf = t2_cleandf[pred]
#     for ec, dist_i in subdf.items():
#         try:
#             j = all_label.index(ec)
#         except ValueError:
#             continue
#         dist_i = infer_confidence_gmm(dist_i, gmm_lst)
#         prelabel[i, j] = dist_i
#     if i%500==0:
#         print('i:',i)

# print('prelabel:',prelabel.shape)
# print('prelabel:',prelabel[0][:10])
# # save the prelabel for clean subset
# outf = '/ibex/user/niuk0a/funcarve/cobra/data/result/split70t2_CLEANt2_subset'
# with open(outf+'_tensor.pkl', 'wb') as f:
#     pickle.dump(prelabel, f)
####################### t2_iterdf
# print('t2_cleandf:',t2_iterdf.shape,'head:\n',t2_iterdf.head(),sep='')
# prelabel_iter = torch.zeros(len(pred_names), len(all_label))

# for i, pred in enumerate(pred_names):
#     subdf = t2_iterdf[pred]
#     for ec, dist_i in subdf.items():
#         try:
#             j = all_label.index(ec)
#         except ValueError:
#             continue
#         dist_i = infer_confidence_gmm(dist_i, gmm_lst)
#         prelabel_iter[i, j] = dist_i
#     if i%100==0:
#         print('i:',i,flush=True)


# prelabel = prelabel_iter
# outf = '/ibex/user/niuk0a/funcarve/cobra/data/result/split70t2_ITERt2_parr'
# with open(outf+'_tensor.pkl', 'wb') as f:
#     pickle.dump(prelabel, f)

# # predlable distance to confidence by gmm
# # Ranking-based metrics
# results = {}
# y_true = label.flatten()
# y_pred_prob = prelabel.flatten()
# results['roc_auc'] = roc_auc_score(y_true, y_pred_prob, average='micro', multi_class='ovr')
# results['average_precision'] = average_precision_score(y_true, y_pred_prob, average='micro')
# results['mse'] = mean_squared_error(y_true, y_pred_prob)
# results['log_loss'] = log_loss(y_true, y_pred_prob)
# results['brier_score'] = brier_score_loss(y_true.flatten(), y_pred_prob.flatten())

# print('results:',results)

## clean result 
# results: {'roc_auc': 0.9929912308631522, 'average_precision': 0.9848405953109239, 'mse': 6.742816e-05, 'log_loss': 0.0006607438641745787, 'brier_score': 6.74281297714513e-05}
## iter result


### evaluates
# results: {'roc_auc': 0.9929913063960483, 'average_precision': 0.9848701036163645, 'mse': 6.741581e-05, 'log_loss': 0.0006600589291368151, 'brier_score': 6.741582316421824e-05}
# results: {'roc_auc': 0.9929913063960483, 'average_precision': 0.9848701036163645, 'mse': 6.7415815e-05, 'log_loss': 0.0006600589291368148, 'brier_score': 6.74158231642183e-05}

##################################################################














# gmm=None
# out_filename ='/ibex/user/niuk0a/funcarve/cobra/data/result/split70t2_CLEANt2_subset'
# # t2_cleandf = t2_cleandf[multicols]
# # write_max_sep_choices(t2_cleandf, out_filename, gmm=gmm)
# pred_names,pred_label = get_pred_labels(out_filename, pred_type='_maxsep')
# print('pred_names:',len(pred_names))
# cleannames = pred_names
# pred_probs = get_pred_probs(out_filename, pred_type='_maxsep')
# print('pred_label:',len(pred_label),'pred_probs:',len(pred_probs))

# true_label,all_label = get_true_labels_mapname(pred_names,'/ibex/user/niuk0a/funcarve/cobra/split70_bactaxo_test2')
# # true_label, all_label = get_true_labels('/ibex/user/niuk0a/funcarve/cobra/split70_bactaxo_test2')

# print('true_label:',len(true_label),'all_label',len(all_label))
# pre, rec, f1, roc, acc = get_eval_metrics(
#     pred_label, pred_probs, true_label, all_label)
# print("############ EC calling results using maximum separation CLEANT2############")
# print('-' * 75)
# print(f'>>> total samples: {len(true_label)} | total ec: {len(all_label)} \n'
#     f'>>> precision: {pre:.5} | recall: {rec:.5}'
#     f'| F1: {f1:.5} | AUC: {roc:.5} ')
# print('-' * 75)
# multicols: 4332
# pred_names: 12408
# pred_label: 12408 pred_probs: 12408
# true_label: 12408 all_label 1362
# ############ EC calling results using maximum separation CLEANT2############
# ---------------------------------------------------------------------------
# >>> total samples: 12408 | total ec: 1362 
# >>> precision: 0.974 | recall: 0.965| F1: 0.966 | AUC: 0.983 
# ---------------------------------------------------------------------------
############ EC calling results using maximum separation CLEANT2############
# ---------------------------------------------------------------------------
# >>> total samples: 4332 | total ec: 696 
# >>> precision: 0.981 | recall: 0.972| F1: 0.974 | AUC: 0.986 
# ---------------------------------------------------------------------------
# exit()

# of = open('/ibex/user/niuk0a/funcarve/cobra/data/result/split70t2_prs_v4.txt','w')

# out_filename = '/ibex/user/niuk0a/funcarve/cobra/data/result/split70t2_ITERt2v4'
# # out_filename = '/ibex/user/niuk0a/funcarve/cobra/data/result/split70t2_ITERt2v4_scale'
# # write_max_sep_choices(t2_iterdf, out_filename, gmm=gmm)

# pred_names,pred_label = get_pred_labels(out_filename, pred_type='_maxsep')
# print('pred_names:',len(pred_names))
# print('intersection:',len([i for i in pred_names if i in cleannames]))  # 0

# pred_probs = get_pred_probs(out_filename, pred_type='_maxsep')
# print('pred_label:',len(pred_label),'pred_probs:',len(pred_probs))

# interprs = [i for i in pred_names if i in cleannames]
# print('intersection:\n',"\n".join(interprs),sep='',file=of)
# cleanprs = [i for i in cleannames if i not in interprs]
# print('cleanprs not in inter:\n',"\n".join(cleanprs),sep='',file=of)
# iterp = [i for i in pred_names if i not in cleannames]
# print('iterp not in clean:\n',"\n".join(iterp),sep='',file=of)

# true_label, all_label = get_true_labels_mapname(pred_names,'/ibex/user/niuk0a/funcarve/cobra/split70_bactaxo_test2')
# of.close()
# pre, rec, f1, roc, acc = get_eval_metrics(
#     pred_label, pred_probs, true_label, all_label)
# print("############ EC calling results using maximum separation ITERT2############")
# print('-' * 75)
# print(f'>>> total samples: {len(true_label)} | total ec: {len(all_label)} \n'
#     f'>>> precision: {pre:.5} | recall: {rec:.5}'
#     f'| F1: {f1:.5}| AUC: {roc:.5} ')
# print('-' * 75)
# # multicols: 4332
# # pred_names: 2503
# # pred_names: 4332
# # intersection: 2379
# # pred_label: 4332 pred_probs: 4332
# # ############ EC calling results using maximum separation ITERT2############
# # ---------------------------------------------------------------------------
# # >>> total samples: 4332 | total ec: 696 
# # >>> precision: 0.981 | recall: 0.972| F1: 0.974 | AUC: 0.986 
# # ---------------------------------------------------------------------------