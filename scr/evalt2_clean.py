import pickle
import pandas as pd
import numpy as np
import os
import csv

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
# data,prcount = get_all_data(prs, folder)
# # save data
# with open('/ibex/user/niuk0a/funcarve/cobra/data/tmp/up_data.pkl', 'wb') as f:
#     pickle.dump(data, f)
# with open('/ibex/user/niuk0a/funcarve/cobra/data/tmp/up_prcount.pkl', 'wb') as f:
#     pickle.dump(prcount, f)


# /ibex/user/niuk0a/funcarve/cobra/data/tmp/UP000002156_340099_differ_predscore_2.pkl
# /ibex/user/niuk0a/funcarve/cobra/data/tmp/UP000002156_340099_ORIpredscore_0.pkl
#load data
# with open('/ibex/user/niuk0a/funcarve/cobra/data/tmp/up_data.pkl', 'rb') as f:
#     data = pickle.load(f)
# with open('/ibex/user/niuk0a/funcarve/cobra/data/tmp/up_prcount.pkl', 'rb') as f:
#     prcount = pickle.load(f)
# print('data',len(data.keys()),list(data.keys())[0],len(data[list(data.keys())[0]].keys()))
# print('wait...')


# combine_data = {}

# for index in data.keys():
#     for pr in data[index].keys():
#         if pr in combine_data:
#             for ec in data[index][pr].keys():
#                 combine_data[pr][ec] = combine_data[pr][ec] + data[index][pr][ec]/prcount[pr]    
#         else:
#             combine_data[pr] = data[index][pr]
#             for ec in combine_data[pr].keys():
#                 combine_data[pr][ec] = combine_data[pr][ec]/prcount[pr]
        
# df = pd.DataFrame(combine_data)
# print('df',df.shape,df.head())
# df  = df.T
# ## save data
# with open('/ibex/user/niuk0a/funcarve/cobra/data/tmp/up_combine_data.pkl', 'wb') as f:
#     pickle.dump(df, f)

## load data
# with open('/ibex/user/niuk0a/funcarve/cobra/data/tmp/up_combine_data.pkl', 'rb') as f:
#     df = pickle.load(f)
# print('df',df.shape,df.head())

t2_cleandf = pd.read_pickle('/ibex/user/niuk0a/CLEAN/app/results/inputs/split70_bactaxo_test2_maxsep_df.pkl')
# /ibex/user/niuk0a/CLEAN/app/results/inputs/split70_bactaxo_test2_maxsep_df.pkl
# /ibex/user/niuk0a/CLEAN/app/results/inputs/split70_bactaxo_test2_maxsep.csv
t2_iterdf = pd.read_pickle('/ibex/user/niuk0a/funcarve/cobra/data/tmp/up_combine_data.pkl')
t2_iterdf = t2_iterdf.T
# print('t2_cleandf:',t2_cleandf.shape,'head:\n',t2_cleandf.head(),sep='')
# print('t2_iterdf:',t2_iterdf.shape,'head:\n',t2_iterdf.head(),sep='')



cleancols = t2_cleandf.columns.to_list()
itercols = t2_iterdf.columns.to_list()
    
multicols = [i for i in itercols if i in prs]
print('multicols:',len(multicols))

t2_cleandf = t2_cleandf[multicols]
t2_iterdf = t2_iterdf[multicols]
print('t2_cleandf:',t2_cleandf.shape,'head:\n',t2_cleandf.head(),sep='')
print('t2_iterdf:',t2_iterdf.shape,'head:\n',t2_iterdf.head(),sep='')

#### scale the data

# mask = (t2_iterdf != t2_cleandf)
# print('type:',type(mask),mask.shape)
# # print('mask:',mask)

# t2_iterdf[mask] = t2_iterdf[mask] -0.5
# print ('t2_iterdf:',t2_iterdf.shape,'head:\n',t2_iterdf.head(),sep='')


########################



t2_cleandf = t2_cleandf[sorted(t2_cleandf.columns)]
t2_iterdf = t2_iterdf[sorted(t2_iterdf.columns)]
from eval_utils import get_eval_metrics, get_pred_labels, get_pred_probs, write_max_sep_choices

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

gmm=None
out_filename ='/ibex/user/niuk0a/funcarve/cobra/data/result/split70t2_CLEANt2_subset'
# t2_cleandf = t2_cleandf[multicols]
# write_max_sep_choices(t2_cleandf, out_filename, gmm=gmm)
pred_names,pred_label = get_pred_labels(out_filename, pred_type='_maxsep')
print('pred_names:',len(pred_names))
cleannames = pred_names
pred_probs = get_pred_probs(out_filename, pred_type='_maxsep')
print('pred_label:',len(pred_label),'pred_probs:',len(pred_probs))

true_label,all_label = get_true_labels_mapname(pred_names,'/ibex/user/niuk0a/funcarve/cobra/split70_bactaxo_test2')
# true_label, all_label = get_true_labels('/ibex/user/niuk0a/funcarve/cobra/split70_bactaxo_test2')

print('true_label:',len(true_label),'all_label',len(all_label))
pre, rec, f1, roc, acc = get_eval_metrics(
    pred_label, pred_probs, true_label, all_label)
print("############ EC calling results using maximum separation CLEANT2############")
print('-' * 75)
print(f'>>> total samples: {len(true_label)} | total ec: {len(all_label)} \n'
    f'>>> precision: {pre:.5} | recall: {rec:.5}'
    f'| F1: {f1:.5} | AUC: {roc:.5} ')
print('-' * 75)
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

of = open('/ibex/user/niuk0a/funcarve/cobra/data/result/split70t2_prs_v4.txt','w')

out_filename = '/ibex/user/niuk0a/funcarve/cobra/data/result/split70t2_ITERt2v4'
# out_filename = '/ibex/user/niuk0a/funcarve/cobra/data/result/split70t2_ITERt2v4_scale'
# write_max_sep_choices(t2_iterdf, out_filename, gmm=gmm)

pred_names,pred_label = get_pred_labels(out_filename, pred_type='_maxsep')
print('pred_names:',len(pred_names))
print('intersection:',len([i for i in pred_names if i in cleannames]))  # 0

pred_probs = get_pred_probs(out_filename, pred_type='_maxsep')
print('pred_label:',len(pred_label),'pred_probs:',len(pred_probs))

interprs = [i for i in pred_names if i in cleannames]
print('intersection:\n',"\n".join(interprs),sep='',file=of)
cleanprs = [i for i in cleannames if i not in interprs]
print('cleanprs not in inter:\n',"\n".join(cleanprs),sep='',file=of)
iterp = [i for i in pred_names if i not in cleannames]
print('iterp not in clean:\n',"\n".join(iterp),sep='',file=of)

true_label, all_label = get_true_labels_mapname(pred_names,'/ibex/user/niuk0a/funcarve/cobra/split70_bactaxo_test2')
of.close()
pre, rec, f1, roc, acc = get_eval_metrics(
    pred_label, pred_probs, true_label, all_label)
print("############ EC calling results using maximum separation ITERT2############")
print('-' * 75)
print(f'>>> total samples: {len(true_label)} | total ec: {len(all_label)} \n'
    f'>>> precision: {pre:.5} | recall: {rec:.5}'
    f'| F1: {f1:.5}| AUC: {roc:.5} ')
print('-' * 75)
# multicols: 4332
# pred_names: 2503
# pred_names: 4332
# intersection: 2379
# pred_label: 4332 pred_probs: 4332
# ############ EC calling results using maximum separation ITERT2############
# ---------------------------------------------------------------------------
# >>> total samples: 4332 | total ec: 696 
# >>> precision: 0.981 | recall: 0.972| F1: 0.974 | AUC: 0.986 
# ---------------------------------------------------------------------------