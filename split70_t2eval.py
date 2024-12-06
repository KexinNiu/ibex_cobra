import pandas as pd
import pickle
import glob
import os
import torch
import csv
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.metrics import precision_score, recall_score, \
    roc_auc_score, accuracy_score, f1_score, average_precision_score
import numpy as np

def make_lastiterdict(dr):
    lastiterdict={}
    for file in os.listdir(dr):
        if str(file).startswith('UP'):
            # print(file)
            upid = file.split('_')[0]
            iterid = file.split('_')[-1].split('.')[0]
            if upid not in lastiterdict:
                lastiterdict[upid] = [int(iterid),file]
            else:
                if int(iterid) > lastiterdict[upid][0]:
                    lastiterdict[upid]=[int(iterid),file]
    print('len(lastiterdict)',len(lastiterdict))  #  12408
    # save lastiterdict
    with open('/ibex/user/niuk0a/funcarve/cobra/split70t2_lastiterdict.pkl','wb') as f:
        pickle.dump(lastiterdict,f)
    return lastiterdict

def load_lastiterdict(path):
    """加载 lastiterdict 数据"""
    with open(path, 'rb') as f:
        lastiterdict = pickle.load(f)
    print('Loaded lastiterdict with keys:', len(lastiterdict))
    return lastiterdict

def process_up_file(file_path, prs, evaldf, times,allprs):
    """处理单个 UP 文件，更新 evaldf 和计数器 times"""
    with open(file_path, 'rb') as f:
        up = pickle.load(f)
        updf = pd.DataFrame.from_dict(up)
        cols = set(updf.columns.to_list())
        allprs.update(cols)
        print('allpres:',len(allprs))   
    # 如果 evaldf 为空，初始化索引
    if evaldf.empty:
        evaldf = pd.DataFrame(index=updf.index)
    
    # 用于收集需要插入的列
    new_columns = {}
    
    for col in updf.columns:
        if type(col) != str:
            continue
        try:
            name = col.split('|')[1]  # 提取列名中间的部分
        except IndexError:
            print(f'Invalid column name "{col}" in file "{file_path}"')
            continue
        
        if name in prs:
            # 更新计数
            times[name] += 1
            
            # 收集新列数据
            if name not in evaldf.columns:
                new_columns[name] = updf[col]
            else:
                evaldf[name] += updf[col]
    
    # 将新列一次性加入 evaldf
    if new_columns:
        new_columns_df = pd.DataFrame(new_columns, index=evaldf.index)
        evaldf = pd.concat([evaldf, new_columns_df], axis=1)

    return evaldf, times,allprs

def average_evaldf(evaldf, times):
    """计算每列的平均值"""
    for col in evaldf.columns:
        if times[col] > 0:
            evaldf[col] = evaldf[col] / times[col]
        else:
            print(f"Warning: No data found for column '{col}'")
    return evaldf

def save_evaldf(evaldf, save_path):
    """保存 evaldf 到指定路径"""
    evaldf.to_pickle(save_path)
    print(f'Evaldf saved to {save_path}') 

def gene_evaldf(lastdicf,savepath='/ibex/user/niuk0a/funcarve/cobra/split70t2_evaldf.pkl'):
    evaldf = pd.DataFrame()
    allprs = set()
    times={key:val for key,val in zip(prs,[0]*len(prs))}
    lastiterdict = load_lastiterdict(lastdicf)
    print('len(lastiterdict)',lastiterdict.keys()) #  43
    # 遍历文件目录
    for file in os.listdir(dr):
        if not file.startswith('UP'):
            continue
        upid = file.split('_')[0]
        iterid = file.split('_')[-1].split('.')[0]
        
        # 检查文件是否为最新迭代
        if lastiterdict.get(upid, [None, None])[1] == file:
            file_path = os.path.join(dr, file)
            print(f'Processing file: {file}')
            evaldf, times, allprs= process_up_file(file_path, prs, evaldf, times,allprs)
    # 计算平均值

    evaldf = average_evaldf(evaldf, times)

    # 保存结果
    # save_evaldf(evaldf, '/ibex/user/niuk0a/funcarve/cobra/split70t2_evaldf.pkl')
    # save_evaldf(evaldf, '/ibex/user/niuk0a/funcarve/cobra/split70t2_evaldf_addon.pkl')
    save_evaldf(evaldf, savepath)

    # 打印基本信息
    print(evaldf.head())
    print(f'Shape of evaldf: {evaldf.shape}')
    evalpres = evaldf.columns.to_list()
    print('missing:',len([i for i in prs if i not in evalpres]))  # 0
    print('allprs:',len(allprs))  # 12408   
    return evaldf

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

def infer_confidence_gmm(distance, gmm_lst):
    confidence = []
    for j in range(len(gmm_lst)):
        main_GMM = gmm_lst[j]
        a, b = main_GMM.means_
        true_model_index = 0 if a[0] < b[0] else 1
        certainty = main_GMM.predict_proba([[distance]])[0][true_model_index]
        confidence.append(certainty)
    return np.mean(confidence)

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

def get_true_labels(file_name):
    result = open(file_name+'.csv', 'r')
    csvreader = csv.reader(result, delimiter='\t')
    all_label = set()
    true_label_dict = {}
    header = True
    count = 0
    for row in csvreader:
        # don't read the header
        if header is False:
            count += 1
            true_ec_lst = row[1].split(';')
            true_label_dict[row[0]] = true_ec_lst
            for ec in true_ec_lst:
                all_label.add(ec)
        if header:
            header = False
    true_label = [true_label_dict[i] for i in true_label_dict.keys()]
    return true_label, all_label

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

def get_pred_probs(out_filename, pred_type="_maxsep"):
    file_name = out_filename+pred_type
    result = open(file_name+'.csv', 'r')
    csvreader = csv.reader(result, delimiter=',')
    pred_probs = []
    for row in csvreader:
        preds_ec_lst = []
        preds_with_dist = row[1:]
        probs = torch.zeros(len(preds_with_dist))
        count = 0
        for pred_ec_dist in preds_with_dist:
            # get EC number 3.5.2.6 from EC:3.5.2.6/10.8359
            ec_i = float(pred_ec_dist.split(":")[1].split("/")[1])
            probs[count] = ec_i
            #preds_ec_lst.append(probs)
            count += 1
        # sigmoid of the negative distances 
        probs = (1 - torch.exp(-1/probs)) / (1 + torch.exp(-1/probs))
        probs = probs/torch.sum(probs)
        pred_probs.append(probs)
    return pred_probs

f='/ibex/user/niuk0a/funcarve/cobra/split70_bactaxo_test2.csv'
df = pd.read_csv(f,sep='\t')
prs = df['Entry'].values.tolist()  # 12408
print('prs',prs[:10],len(prs))  # 12408

# multidf = pd.read_csv('/ibex/user/niuk0a/CLEAN/app/data/taxodata/split70_bactaxo_test2_multilabel.csv',sep='\t')
# prs = multidf['Entry'].values.tolist()  # 627

dr = '/ibex/user/niuk0a/funcarve/cobra/'
# lastiterdict = make_lastiterdict(dr)
# evaldf = pd.DataFrame()

# evaldf = gene_evaldf('/ibex/user/niuk0a/funcarve/cobra/split70t2_lastiterdict.pkl','/ibex/user/niuk0a/funcarve/cobra/split70t2_evaldf_addon.pkl')

## load evaldf
# with open('/ibex/user/niuk0a/funcarve/cobra/split70t2_evaldf.pkl','rb') as f:
#     evaldf = pickle.load(f)
# with open('/ibex/user/niuk0a/funcarve/cobra/split70t2_evaldf_addon.pkl','rb') as f:
#     evaldf = pickle.load(f)
# print('evaldf:',evaldf.shape)
# evalpres = evaldf.columns.to_list()
# print('missing:',len([i for i in prs if i not in evalpres]))  # 2414/12408


t2_cleandf = pd.read_pickle('/ibex/user/niuk0a/CLEAN/app/results/inputs/split70_bactaxo_test2_maxsep_df.pkl')
# t2_iterdf = pd.read_pickle('/ibex/user/niuk0a/funcarve/cobra/split70t2_evaldf.pkl')
t2_iterdf = pd.read_pickle('/ibex/user/niuk0a/funcarve/cobra/split70t2_evaldf_addon.pkl')
print('t2_cleandf:',t2_cleandf.shape)
# print('t2_iterdf:',t2_iterdf.shape)
print('t2addon_iterdf:',t2_iterdf.shape)

# OREDER THE COLUMNS 
# t2_cleandf = t2_cleandf[sorted(t2_cleandf.columns)]
# t2_iterdf = t2_iterdf[sorted(t2_iterdf.columns)]

print('t2_cleandf:',t2_cleandf.shape,'head:\n',t2_cleandf.head(),sep='')
print('t2_iterdf:',t2_iterdf.shape,'head:\n',t2_iterdf.head(),sep='')

cleancols = t2_cleandf.columns.to_list()
itercols = t2_iterdf.columns.to_list()

multicols = [i for i in itercols if i in prs]
print('multicols:',len(multicols))

t2_cleandf = t2_cleandf[multicols]
t2_iterdf = t2_iterdf[multicols]
# smallt2_cleandf = t2_cleandf[itercols]
# t2_cleandf = smallt2_cleandf
# print('t2_cleandf:',t2_cleandf.shape)

# sort columns
t2_cleandf = t2_cleandf[sorted(t2_cleandf.columns)]
t2_iterdf = t2_iterdf[sorted(t2_iterdf.columns)]
# for i in range(len(t2_cleandf.columns)):
#     name = t2_cleandf.columns[i]


# print(t2_cleandf[name].head())
print('t2_iterdf',t2_iterdf.shape)
        


gmm=None
out_filename = '/ibex/user/niuk0a/funcarve/cobra/split70t2_CLEANt2_multi'
# write_max_sep_choices(t2_cleandf, out_filename, gmm=gmm)
pred_names,pred_label = get_pred_labels(out_filename, pred_type='_maxsep')
print('pred_names:',len(pred_names))
cleannames = pred_names
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
#     f'>>> precision: {pre:.3} | recall: {rec:.3}'
#     f'| F1: {f1:.3} | AUC: {roc:.3} ')
# print('-' * 75)

############ EC calling results using maximum separation CLEANT2############
# ---------------------------------------------------------------------------
# >>> total samples: 627 | total ec: 227 
# >>> precision: 0.934 | recall: 0.836| F1: 0.857 | AUC: 0.918 
# ---------------------------------------------------------------------------
############ EC calling results using maximum separation CLEANT2############
# ---------------------------------------------------------------------------
# >>> total samples: 2503 | total ec: 548 
# >>> precision: 0.986 | recall: 0.978| F1: 0.98 | AUC: 0.989 
# ---------------------------------------------------------------------------
############ EC calling results using maximum separation CLEANT2############
# ---------------------------------------------------------------------------
# >>> total samples: 142 | total ec: 75 
# >>> precision: 0.962 | recall: 0.896| F1: 0.911 | AUC: 0.948 
# ---------------------------------------------------------------------------


of = open('/ibex/user/niuk0a/funcarve/cobra/split70t2_prs_v3.txt','w')

out_filename = '/ibex/user/niuk0a/funcarve/cobra/split70t2_ITERt2v3'
write_max_sep_choices(t2_iterdf, out_filename, gmm=gmm)

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
    f'>>> precision: {pre:.3} | recall: {rec:.3}'
    f'| F1: {f1:.3} | AUC: {roc:.3} ')
print('-' * 75)
exit()
############ EC calling results using maximum separation ITERT2############
# ---------------------------------------------------------------------------
# >>> total samples: 613 | total ec: 281 
# >>> precision: 0.993 | recall: 0.991| F1: 0.99 | AUC: 0.995 
# ---------------------------------------------------------------------------

# of = open('/ibex/user/niuk0a/funcarve/cobra/split70t2_prs_addon_multi.txt','w')

# out_filename = '/ibex/user/niuk0a/funcarve/cobra/split70t2_ITERt2addon_multi'
# write_max_sep_choices(t2_iterdf, out_filename, gmm=gmm)

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
# print("############ EC calling results using maximum separation ITERT2 !addon ############")
# print('-' * 75)
# print(f'>>> total samples: {len(true_label)} | total ec: {len(all_label)} \n'
#     f'>>> precision: {pre:.3} | recall: {rec:.3}'
#     f'| F1: {f1:.3} | AUC: {roc:.3} ')
# print('-' * 75)

# ############ EC calling results using maximum separation ITERT2 !addon ############
# ---------------------------------------------------------------------------
# >>> total samples: 2503 | total ec: 548 
# >>> precision: 0.986 | recall: 0.978| F1: 0.98 | AUC: 0.989 
############ EC calling results using maximum separation ITERT2 !addon multi############
# ---------------------------------------------------------------------------
# >>> total samples: 142 | total ec: 75 
# >>> precision: 0.962 | recall: 0.896| F1: 0.911 | AUC: 0.948 
# ---------------------------------------------------------------------------