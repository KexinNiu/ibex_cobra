import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# f='/ibex/user/niuk0a/funcarve/cobra/eval_clean_all_7_allec.csv'

#############################################
################ec uniprot###############
#############################################

# 读取数据
# file_path = '/ibex/user/niuk0a/funcarve/cobra/eval_7results.csv'  # 替换为你的文件路径
# file_path ='/ibex/user/niuk0a/funcarve/cobra/eval_clean_all_7result_add.csv'
file_path ='/ibex/user/niuk0a/funcarve/cobra/eval_clean_all_7result_add1.csv'
# file_path ='/ibex/user/niuk0a/funcarve/cobra/eval_clean_7v3.csv'

data = pd.read_csv(file_path)

print (data.head())
## for each type with different iter_round, we only keeps the one with the highest roc
data = data.sort_values(by='roc', ascending=False).drop_duplicates(subset=['modelname', 'type','thr'], keep='first')

print (data.head())

# 确保数据类型正确
# data['thr'] = data['thr'].astype(float)
data['roc'] = data['roc'].astype(float)
# reorder by modelname then type and thr
data = data.sort_values(by=['modelname', 'type', 'thr'])

# 获取所有的模型名称
models = data['modelname'].unique()

# 设置大图网格
n_rows = (len(models) + 3) // 4  # 每行最多放4个小图
fig, axes = plt.subplots(n_rows, 4, figsize=(20, 5 * n_rows), squeeze=False)

# 遍历每个模型，生成对应的小图
for idx, model in enumerate(models):
    ax = axes[idx // 4][idx % 4]  # 选择子图
    model_data = data[data['modelname'] == model]
    # keep thres < 1
    # model_data = model_data[model_data['thr'] < 1]

    # 绘制每种 type 的曲线
    for type_name, type_data in model_data.groupby('type'):
        if type_name == 'CLEAN':
            continue
        ax.plot(type_data['thr'], type_data['roc'], label=type_name, marker='o')

    # 添加水平参考线
    clean_roc = model_data[model_data['type'] == 'CLEAN']['roc'].values
    print('clean_roc=',clean_roc)
    if clean_roc.size > 0:
        clean_value = clean_roc[0]
        ax.axhline(y=clean_value, color='gray', linestyle='--', label='CLEAN Reference')
    
    # tiaozhaeng y轴范围
    # ax.set_xlim(0, 5)
    ax.set_ylim(clean_value-0.01, model_data['roc'].max() + 0.01)

    # 设置图例和标题
    ax.set_title(f'Model: {model}')
    ax.set_xlabel('Penalty score')
    ax.set_ylabel('ROC')
    # ax.legend()

# 移除多余的空白子图
for j in range(len(models), n_rows * 4):
    fig.delaxes(axes[j // 4][j % 4])

# 调整布局并保存
# add legend for all subplots by using the last one and put it outside,at right bottom
plt.legend(loc='best', bbox_to_anchor=(1, 1), ncol=1)
#'best', 'upper right', 'upper left', 'lower left', 'lower right', 'right', 'center left', 'center right', 'lower center', 'upper center', 'center'

plt.tight_layout()
# plt.savefig('eval_thr_plotv3.png', dpi=300)
plt.savefig('eval_thr_plot.png', dpi=300)

plt.show()

exit()

######################################################################################################################3
#############################################
################model coverage###############
#############################################
# 读取数据
# file_path = '/ibex/user/niuk0a/funcarve/cobra/model_bigg_overlap_i1.csv'  # 替换为你的文件路径
# data = pd.read_csv(file_path)
# ## data add a new column precision
# data['precision'] = data['both']/(data['predmodel'])
# data['recall'] = data['both']/(data['biggmodel'])
# data['f1'] = 2*data['precision']*data['recall']/(data['precision']+data['recall'])

# file_path ='/ibex/user/niuk0a/funcarve/cobra/emodel_bigg_overlap_rerun.csv'
file_path ='/ibex/user/niuk0a/funcarve/cobra/emodel_bigg_overlap_v3.csv'
df = pd.read_csv(file_path)

recon_file_path = '/ibex/user/niuk0a/funcarve/cobra/model_bigg_overlap_i1_recon.csv'
recon_data = pd.read_csv(recon_file_path)

# data = pd.concat([df, data], ignore_index=True)
data = df
# remove thr =  0
data = data[data['thr'] != '0']
data = data[data['iter'] < 5]


data['title'] = data['modelname']
# 	genome 	ORGID	ORGNASIM
# iLJ478	NC_000853.1	243274	Thermotoga maritima MSB8
# iJN1463	NC_002947.4	160488	Pseudomonas putida KT2440
# iCN900	AM180355.1	272563	Clostridioides difficile 630
# iAF1260	NC_000913.3	83333	Escherichia coli str. K-12 substr. MG1655
# iAF987	CP000148.1	269799	Geobacter metallireducens GS-15
# iNJ661	NC_000962.3	83332	Mycobacterium tuberculosis H37Rv
# iYO844	AL009126.3	224308	Bacillus subtilis subsp. subtilis str. 168
# rename the title to the organism name
for index, row in data.iterrows():
    if row['modelname'] == 'iLJ478':
        data.at[index, 'title'] = 'Thermotoga maritima MSB8'
    if row['modelname'] == 'iJN1463':
        data.at[index, 'title'] = 'Pseudomonas putida KT2440'
    if row['modelname'] == 'iCN900':
        data.at[index, 'title'] = 'Clostridioides difficile 630'
    if row['modelname'] == 'iAF1260':
        data.at[index, 'title'] = 'Escherichia coli str. K-12 substr. MG1655'
    if row['modelname'] == 'iAF987':
        data.at[index, 'title'] = 'Geobacter metallireducens GS-15'
    if row['modelname'] == 'iNJ661':
        data.at[index, 'title'] = 'Mycobacterium tuberculosis H37Rv'
    if row['modelname'] == 'iYO844':
        data.at[index, 'title'] = 'Bacillus subtilis subsp. subtilis str. 168'


# print (data.head())
## for each type with different iter_round, we only keeps the one with the highest roc
# data = data.sort_values(by='roc', ascending=False).drop_duplicates(subset=['modelname', 'evaltype','thr','iter'], keep='first')

# evalcol = 'f1'
# evalcol='tp'
evalcol='ratio_biggmodel'    

# 设置绘图风格
sns.set(style="whitegrid")

# 创建 FacetGrid，每个 modelname 一个小图
g = sns.FacetGrid(
    data,
    # col="modelname",  # 每个 modelname 一列小图
    col="title",  # 每个 modelname 一列小图

    col_wrap=4,       # 每行放 4 个小图，可以根据需要调整
    height=4,         # 每个小图的高度
    sharey=False       # 所有小图共享 y 轴
)

# 添加线条图到每个子图中
g.map_dataframe(
    sns.lineplot, 
    x="iter", 
    # y="both", 
    # y="ratio_biggmodel", 
    # y="ratio_predmodel", 
    y=evalcol, 
    hue="thr",  # 用 thr 区分线条
    style="thr", 
    markers=True, 
    dashes=False
)
#给每个字图加一条水平参考线 数据来自recon_data
# g.map(plt.axhline, y=recon_data[evalcol].values[0], ls='--', c='gray')


for ax, modelname in zip(g.axes.flat, data["modelname"].unique()):
    #添加参考线
    # print('val=',recon_data[recon_data['modelname']==modelname][evalcol].values[0])
    refval = recon_data[recon_data['modelname']==modelname][evalcol].values[0]
    # ax.axhline(y=refval, color='grey', linestyle='--', label='Recon {refval}')
    #set legend for recon reference
    # ax.legend([f'Recon={refval}'], loc='lower left',color='grey')
    ax.text(1, 0.05, f'Recon={refval:.3f}', horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes, color='black', fontsize=10)
    subset = data[data["modelname"] == modelname]  # 提取当前子图的对应数据
    # ymin, ymax = subset["both"].min(), subset["both"].max()  # 获取当前数据的 y 值范围
    # ymin, ymax = subset["ratio_biggmodel"].min(), subset["ratio_biggmodel"].max()  # 获取当前数据的 y 值范围
    ymin, ymax = subset["ratio_predmodel"].min(), subset["ratio_predmodel"].max()  # 获取当前数据的 y 值范围
    ymin, ymax = subset[evalcol].min(), subset[evalcol].max()  # 获取当前数据的 y 值范围
    # ax.set_ylim(ymin - 5, ymax + 5)  # 设置 y 轴范围，增加少量缓冲区以美观
    # ax.set_ylim(ymin-0.01, ymax + 0.02)  # 设置 y 轴范围，增加少量缓冲区以美观
    ax.set_ylim(ymin-0.01, ymax + 0.01)  # 设置 y 轴范围，增加少量缓冲区以美观
    # ax.set_ylim(refval-0.001, ymax + 0.01)  # 设置 y 轴范围，增加少量缓冲区以美观
    # ax.set_ylim(refval-1, ymax + 50)  # 设置 y 轴范围，增加少量缓冲区以美观
# 设置图例和标题
g.add_legend(title="Threshold")
## set legend for recon reference
# g.add_legend(title="Threshold", labels=['Recon Reference'], loc='best')
# g.set_titles("{col_name}")  # 每个子图的标题为对应的 modelname
# 字图的标题为dataframe中的title
# print(data['title'].unique())
# print('col_name=',g.col_names)
g.set_titles("{col_name}")  # 每个子图的标题为对应的 modelname
# g.set_axis_labels("Iteration", "Both")
# g.set_axis_labels("Iteration", "ratio_biggmodel")
# g.set_axis_labels("Iteration", "ratio_predmodel")
g.set_axis_labels("Iteration", evalcol)

# 调整整体布局
plt.subplots_adjust(top=0.9)
# g.fig.suptitle("Both vs Iter for Each Modelname", fontsize=16)  # 总标题
# g.fig.suptitle("ratio_biggmodel vs Iter for Each Modelname", fontsize=16)  # 总标题
# g.fig.suptitle("ratio_predmodel vs Iter for Each Modelname", fontsize=16)  # 总标题
g.fig.suptitle(f"Reactions {evalcol} based on BiGG models ", fontsize=16)  # 总标题

# 显示图表
plt.show()
# plt.savefig('eval_modelrxn_ratio_biggmodel.png', dpi=300)
# plt.savefig('eval_modelrxn_ratio_predmodel.png', dpi=300)
# plt.savefig(f'eval_modelrxn_{evalcol}_recon5.png', dpi=300)
plt.savefig(f'eval_modelrxn_{evalcol}_v3.png', dpi=300)

exit()


#############################################
################model memote  ###############
#############################################

# def eval_memote()