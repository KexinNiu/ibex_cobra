from eval_utils import *
import pandas as pd
import os


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
resultmatrix = {
    'modelname':[],
    'type':[],
    'iter_round':[],
    'thr':[],
    'reward':[],    
    'pre':[],
    'rec':[],
    'f1':[],
    'roc':[],
    'acc':[],
    'pred_len':[],
    'all_len':[]
}
evallist=['allec','flux','fluxblock','block']
rewardlist = [0.2,0.5,1.0,3.0]
threlist = [5,6]

# evallist=['allec']
# rewardlist = [0.2]
# threlist = [5]

def eval_ori(pref,unif):
    pred_label, pred_names, pred_all_label = get_pred_labels(pref)
    pred_probs = get_pred_probs(pref)
    # print('pred names:',len(pred_names),pred_names[:5])
    true_label, all_label, pred_names, pred_label, pred_probs = get_true_labels(unif,pref,pred_names,pred_label,pred_probs,pred_all_label)
    return true_label, all_label, pred_names, pred_label, pred_probs




############################GENERATE EVALUATION RESULT#############################
# data_root='/ibex/user/niuk0a/funcarve/cobra/data/tmp'
# # print('only for iYO844')
# # prefs = prefs[-1:]
# # unifs = unifs[-1:]
# # modelnames = modelnames[-1:]
# # print('prefs:',prefs)
# # print('unifs:',unifs)
# # print('modelnames:',modelnames)



# for i in range(0,len(prefs)):
#     pref = prefs[i]
#     unif = unifs[i]
#     modelname = modelnames[i]
#     print('pref:',pref)
#     print('unif:',unif)
#     print('modelname:',modelname)
#     ### load original CLEAN prediction and dfs
#     true_label, all_label, pred_names, pred_label, pred_probs = eval_ori(pref,unif)
#     ### get original evaluation results
#     pre,rec,f1,roc,acc = get_eval_metrics(pred_label, pred_probs, true_label, all_label)
    
#     resultmatrix['modelname'].append(modelname)
#     resultmatrix['type'].append('CLEAN')
#     resultmatrix['iter_round'].append(0)
#     resultmatrix['thr'].append(0)
#     resultmatrix['reward'].append(0)    
#     resultmatrix['pre'].append(pre)
#     resultmatrix['rec'].append(rec)
#     resultmatrix['f1'].append(f1)
#     resultmatrix['roc'].append(roc)
#     resultmatrix['acc'].append(acc)
#     resultmatrix['pred_len'].append(len(pred_names))
#     resultmatrix['all_len'].append(len(all_label))
#     print(f'{modelname} CLEAN 0 0 0 {pre} {roc} {acc}',flush=True)
#     ### done for original CLEAN prediction

#     # load and evaluate enzyme annotation results for each condition
#     for reward in rewardlist:
#         if reward < 1.0:
#             rewardname = 'p'+str(reward).split('.')[0]+str(reward).split('.')[1]
#         else:
#             rewaardname = str(reward).split('.')[0]

#         for thre in threlist:
#             thrname = str(thre)
#             for evaltype in evallist:
#                 # /ibex/user/niuk0a/funcarve/cobra/data/tmp/iYO844_allec_R3_T5_t5.0_ORIpredscore_0.pkl
#                 orif=f'{data_root}/{modelname}_{evaltype}_R{rewardname}_T{thrname}_ORIpredscore_0.pkl'
#                 # FleNotFoundError: [Errno 2] No such file or directory: '/ibex/user/niuk0a/funcarve/cobra/data/tmp/iLJ478_allec_Rp02_T5_t0.2_ORIpredscore_0.pkl'

#                 orid = pd.read_pickle(orif)
#                 for i in range(1,4):
#                     differf = f'{data_root}/{modelname}_{evaltype}_R{rewardname}_T{thrname}_differ_predscore_{i}.pkl'
#                     if not os.path.exists(differf):
#                         continue
#                     # load and evaluate enzyme annotation results
#                     differd = pd.read_pickle(differf)
#                     currentd = update_differ(orid,differd)
#                     interpred_names,interpred_label_cutoff,interpred_label,interpred_probs = eval_df(currentd,all_label,pred_names,thre)
#                     pre,rec,f1,roc,acc = get_eval_metrics(interpred_label, interpred_probs, true_label, all_label)
#                     resultmatrix['modelname'].append(modelname)
#                     resultmatrix['type'].append(evaltype)
#                     resultmatrix['iter_round'].append(i)
#                     resultmatrix['thr'].append(thre)
#                     resultmatrix['reward'].append(reward)
#                     resultmatrix['pre'].append(pre)
#                     resultmatrix['rec'].append(rec)
#                     resultmatrix['f1'].append(f1)
#                     resultmatrix['roc'].append(roc)
#                     resultmatrix['acc'].append(acc)
#                     resultmatrix['pred_len'].append(len(interpred_names))
#                     resultmatrix['all_len'].append(len(all_label))
#                     print(f'{modelname} {evaltype} {reward} {thre} {i} done',flush=True)


# file_path = '/ibex/user/niuk0a/funcarve/cobra/data/result/eval_clean_7v4_specific_ec.csv'

# df = pd.DataFrame(resultmatrix)
# df.to_csv(file_path,index=False)
# print(f'Evaluation results saved to {file_path}')
# exit()
###############################     PLOTS      ###################################
from matplotlib import pyplot as plt

# file_path = '/ibex/user/niuk0a/funcarve/cobra/data/result/eval_clean_7v4.csv'
file_path = '/ibex/user/niuk0a/funcarve/cobra/data/result/eval_clean_7v4_specific_ec.csv'

data = pd.read_csv(file_path,sep=',')
# data1 = data[(data['type']=='allec')]
# data2 = data[(data['type']=='CLEAN')]
# data = pd.concat([data1,data2])
print (data.head())
# data = data.sort_values(by='roc', ascending=False).drop_duplicates(subset=['modelname', 'type','thr','reward'], keep='first')
data = data.sort_values(by='f1', ascending=False).drop_duplicates(subset=['modelname', 'type','thr','reward'], keep='first')
print (data.head())
modelnames = ['iLJ478','iAF1260','iJN1463','iYO844','iNJ661','iAF987','iCN900']
# taxonames = ['Thermotoga maritima MSB8','Escherichia coli str. K-12 substr. MG1655','Pseudomonas putida KT2440','Bacillus subtilis subsp. subtilis str. 168','Mycobacterium tuberculosis H37Rv','Geobacter metallireducens GS-15','Clostridioides difficile 630']
taxonames = ['T. maritima MSB8','E. coli MG1655','P. putida KT2440','B. subtilis 168','M. tuberculosis H37Rv','G. metallireducens GS-15','C. difficile 630']
# taxonames = ['T. mari','E. coli','P. putida','B. sub','M. tuber','G. meta','C. diff']
# 确保数据类型正确
# data['thr'] = data['thr'].astype(float)
data['roc'] = data['roc'].astype(float)
data['f1'] = data['f1'].astype(float)
data = data.sort_values(by=['modelname', 'type', 'thr','reward'])
evaltype='roc'
models = data['modelname'].unique()

# 设置大图网格
n_rows = (len(models) + 3) // 4  # 每行最多放4个小图
fig, axes = plt.subplots(n_rows, 4, figsize=(20, 5 * n_rows), squeeze=False)

# 遍历每个模型，生成对应的小图
for idx, model in enumerate(models):
    ax = axes[idx // 4][idx % 4]  # 选择子图
    model_data = data[data['modelname'] == model]
    print('model:',model,'len:',len(model_data))
    # 绘制每种 type 的曲线
    for type_name, type_data in model_data.groupby('type'):
        # print('type:',type_data)
        if type_name == 'CLEAN':
            continue
        for thr in type_data['thr'].unique():
            # Filter rows corresponding to the current 'thr'
            thr_data = type_data[type_data['thr'] == thr]
            
            # Set the marker based on 'thr'
            marker = '+' if thr == 5 else 'x'

            # Plot the data with the selected marker
            # ax.plot(thr_data['reward'], thr_data['roc'], 
            ax.plot(thr_data['reward'], thr_data[evaltype], 
                label=f"{type_name} (thr={thr})", marker=marker,alpha=0.7)

    # 添加水平参考线
    # clean_roc = model_data[model_data['type'] == 'CLEAN']['roc'].values
    clean_roc = model_data[model_data['type'] == 'CLEAN'][evaltype].values
    # print('clean_roc=',clean_roc)
    if clean_roc.size > 0:
        clean_value = clean_roc[0]
        ax.axhline(y=clean_value, color='gray', linestyle='--', label='CLEAN Reference')
    
    # tiaozhaeng y轴范围
    # ax.set_xlim(0, 5)
    # ax.set_ylim(clean_value-0.01, model_data['roc'].max() + 0.01)
    ax.set_ylim(min(model_data[evaltype].min() - 0.01,clean_value-0.01), model_data[evaltype].max() + 0.01)

    # 设置图例和标题
    # ax.set_title(f'Model: {model}')
    ax.set_title(f'{taxonames[modelnames.index(model)]}')
    ax.set_xlabel('Rewards score')
    # ax.set_ylabel('ROC')
    if evaltype=='roc':
        ax.set_ylabel('ROC')
    else:
        ax.set_ylabel(evaltype)
    # ax.legend()

# 移除多余的空白子图
for j in range(len(models), n_rows * 4):
    fig.delaxes(axes[j // 4][j % 4])

# 调整布局并保存
# add legend for all subplots by using the last one and put it outside,at right bottom
plt.legend(loc='best', bbox_to_anchor=(1, 1), ncol=1)
#'best', 'upper right', 'upper left', 'lower left', 'lower right', 'right', 'center left', 'center right', 'lower center', 'upper center', 'center'

plt.tight_layout()
# print('save to ../data/result/figures/eval7_enzymeCLEAN_f1_spec.png')
plt.savefig(f'../data/result/figures/eval7_enzymeCLEAN_{evaltype}_spec.png', dpi=300)

# plt.savefig('../data/result/figures/eval7_enzymeCLEAN_roc_spec.png', dpi=300)
plt.show()
exit()
################3plot with increasing iteration####################
# from matplotlib import pyplot as plt
# f = '/ibex/user/niuk0a/funcarve/cobra/data/result/eval_clean_7v4.csv'
# import argparse
# parser = argparse.ArgumentParser(description='Plot evaluation results')
# parser.add_argument('--metric', type=str, default='f1', help='Metric to plot')
# parser.add_argument('--reward', type=float, default=1.0, help='Reward score to plot')
# parser.add_argument('--threshold', type=int, default=5, help='Threshold to plot')
# args = parser.parse_args()

# metric = args.metric
# take_reward = args.reward
# take_thr = args.threshold

# data = pd.read_csv(f,sep=',')


# recondata = data[data['type']=='CLEAN']
# data = data[data['reward'] == take_reward]
# data = data[data['thr'] == take_thr]
# data = pd.concat([recondata,data])
# # data = data.sort_values(by='roc', ascending=False).drop_duplicates(subset=['modelname', 'type','thr','reward'], keep='first')
# # data = data.sort_values(by=metric, ascending=False).drop_duplicates(subset=['modelname', 'type','thr','reward'], keep='first')
# # 确保数据类型正确
# # data['thr'] = data['thr'].astype(float)
# data['roc'] = data['roc'].astype(float)
# data[metric] = data['f1'].astype(float)

# data = data.sort_values(by=['modelname', 'type', 'thr','reward','iter_round'])
# models = data['modelname'].unique()
# # 设置大图网格
# n_rows = (len(models) + 3) // 4  # 每行最多放4个小图
# fig, axes = plt.subplots(n_rows, 4, figsize=(20, 5 * n_rows), squeeze=False)

# # 遍历每个模型，生成对应的小图
# for idx, model in enumerate(models):
#     ax = axes[idx // 4][idx % 4]  # 选择子图
#     model_data = data[data['modelname'] == model]
#     print('model:',model,'len:',len(model_data))
#     # 绘制每种 type 的曲线
#     for type_name, type_data in model_data.groupby('type'):
#         # print('type:',type_data)
#         if type_name == 'CLEAN':
#             continue
#         for thr in type_data['thr'].unique():
#             # Filter rows corresponding to the current 'thr'
#             thr_data = type_data[type_data['thr'] == thr]
            
#             # Set the marker based on 'thr'
#             marker = 'o' if thr == 5 else 'x'

#             # Plot the data with the selected marker
#             # ax.plot(thr_data['reward'], thr_data['roc'], 
#             ax.plot(thr_data['iter_round'], thr_data[metric], 
#                 label=f"{type_name} (thr={thr})", marker=marker)

#     # 添加水平参考线
#     # clean_roc = model_data[model_data['type'] == 'CLEAN']['roc'].values
#     clean_roc = model_data[model_data['type'] == 'CLEAN'][metric].values
#     # print('clean_roc=',clean_roc)
#     if clean_roc.size > 0:
#         clean_value = clean_roc[0]
#         ax.axhline(y=clean_value, color='gray', linestyle='--', label='CLEAN Reference')
    
#     # tiaozhaeng y轴范围
#     # ax.set_xlim(0, 5)
#     # ax.set_ylim(clean_value-0.01, model_data['roc'].max() + 0.01)
#     ax.set_ylim(min(model_data[metric].min() - 0.01,clean_value-0.01), model_data[metric].max() + 0.01)
#     ##shezhi x轴
#     plt.xticks(range(1,3))
#     # 设置图例和标题
#     ax.set_title(f'Model: {model}')
#     ax.set_xlabel('Iteration')
#     # ax.set_ylabel('ROC')
#     ax.set_ylabel(metric)
#     # ax.legend()

# # 移除多余的空白子图
# for j in range(len(models), n_rows * 4):
#     fig.delaxes(axes[j // 4][j % 4])

# # 调整布局并保存
# # add legend for all subplots by using the last one and put it outside,at right bottom
# plt.legend(loc='best', bbox_to_anchor=(1, 1), ncol=1)
# #'best', 'upper right', 'upper left', 'lower left', 'lower right', 'right', 'center left', 'center right', 'lower center', 'upper center', 'center'
# ## add the big title
# fig.suptitle(f'Evaluation of enzyme annotation with reward {take_reward} and threshold {take_thr}', fontsize=16)
# plt.tight_layout()
# plt.savefig(f'../data/result/eval7_enzymeCLEAN_T{take_thr}_R{take_reward}_{metric}.png', dpi=300)
# plt.show()
