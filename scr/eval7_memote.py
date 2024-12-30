from matplotlib import pyplot as plt
import pandas as pd
#################################################

file_path = '/ibex/user/niuk0a/funcarve/cobra/data/result/eval_memote_7v4.csv'

data = pd.read_csv(file_path,sep=',')
# data1 = data[(data['type']=='allec')]
# data2 = data[(data['type']=='CLEAN')]
# data = pd.concat([data1,data2])
print (data.head())
# modelname,evaltype,reward,threshold,iteration,total_metabolites,total_reactions,total_genes,total_compartments,metabolic_coverage,unconserved_metabolites,consistency,annotation_met,annotation_rxn,annotation_gene,annotation_sbo,total_score,newtotal,con_stoi,con_mass,con_charge,con_met,con_unbound,con_unbound_n


modelnames = ['iLJ478','iAF1260','iJN1463','iYO844','iNJ661','iAF987','iCN900']
# taxonames = ['Thermotoga maritima MSB8','Escherichia coli str. K-12 substr. MG1655','Pseudomonas putida KT2440','Bacillus subtilis subsp. subtilis str. 168','Mycobacterium tuberculosis H37Rv','Geobacter metallireducens GS-15','Clostridioides difficile 630']
taxonames = ['T. maritima MSB8','E. coli MG1655','P. putida KT2440','B. subtilis 168','M. tuberculosis H37Rv','G. metallireducens GS-15','C. difficile 630']
# taxonames = ['T. mari','E. coli','P. putida','B. sub','M. tuber','G. meta','C. diff']
# 确保数据类型正确
# only keep evaltype = 'allec'
# data=data[data['evaltype'] == 'allec']

# data['thr'] = data['thr'].astype(float)
data['newtotal'] = data['newtotal'].astype(float)
data['total_score'] = data['total_score'].astype(float)
data = data.sort_values(by=['modelname', 'evaltype', 'threshold','reward','iteration'])
data = data[data['reward'] == '1']
eval='total_score'
eval='total_reactions'
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
    for type_name, type_data in model_data.groupby('evaltype'):
        # print('type:',type_data)
        # if type_name == 'CLEAN':
        #     continue
        # type_data = type_data[type_data['evaltype'] == 'allec']
        type_data = type_data.sort_values(by='iteration')
        # print('type_data:',type_data)

        for thr in type_data['threshold'].unique():
            # Filter rows corresponding to the current 'thr'
            thr_data = type_data[type_data['threshold'] == thr]
            
            # Set the marker based on 'thr'
            # marker = '.' if thr == 5 else 'x'
            marker ='o'
            # Plot the data with the selected marker
            # ax.plot(thr_data['reward'], thr_data['roc'], 
            ax.plot(thr_data['iteration'], thr_data[eval], 
                label=f"{type_name} (threshold={thr})", marker=marker,alpha=0.7)

    # 添加水平参考线
    
    # clean_roc = model_data[model_data['type'] == 'CLEAN'][eval].values

    # if clean_roc.size > 0:
    #     clean_value = clean_roc[0]
    #     ax.axhline(y=clean_value, color='gray', linestyle='--', label='CLEAN Reference')
    

    # tiaozhaeng y轴范围
    # ax.set_xlim(0, 5)
    # ax.set_ylim(clean_value-0.01, model_data['roc'].max() + 0.01)
    # ax.set_ylim(min(model_data[eval].min() - 0.01,clean_value-0.01), model_data[evaltype].max() + 0.01)
    # ax.set_ylim(model_data[eval].min() - 0.01,model_data[eval].max() + 0.01)
    ax.set_ylim(model_data[eval].min() -500,model_data[eval].max() +500)
    # ax.set_ylim(0.79,0.90)

    # 设置图例和标题
    # ax.set_title(f'Model: {model}')
    ax.set_title(f'{taxonames[modelnames.index(model)]}')
    ax.set_xlabel('Iteration')
    # ax.set_ylabel('ROC')
    if eval=='roc':
        ax.set_ylabel('ROC')
    else:
        ax.set_ylabel(eval)
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
plt.savefig(f'../data/result/figures/eval7_memote_{eval}.png', dpi=300)

# plt.savefig('../data/result/figures/eval7_enzymeCLEAN_roc_spec.png', dpi=300)
plt.show()
exit()