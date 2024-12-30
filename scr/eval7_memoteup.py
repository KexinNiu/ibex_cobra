from matplotlib import pyplot as plt
import pandas as pd
#################################################

# file_path = '/ibex/user/niuk0a/funcarve/cobra/data/result/eval_memote_7v4.csv'
file_path = '/ibex/user/niuk0a/funcarve/cobra/data/result/eval_memote_t2up.csv'

data = pd.read_csv(file_path,sep=',')

print (data.head())
# modelname,evaltype,reward,threshold,iteration,total_metabolites,total_reactions,total_genes,total_compartments,metabolic_coverage,unconserved_metabolites,consistency,annotation_met,annotation_rxn,annotation_gene,annotation_sbo,total_score,newtotal,con_stoi,con_mass,con_charge,con_met,con_unbound,con_unbound_n


modelnames = ['iLJ478','iAF1260','iJN1463','iYO844','iNJ661','iAF987','iCN900']
# taxonames = ['Thermotoga maritima MSB8','Escherichia coli str. K-12 substr. MG1655','Pseudomonas putida KT2440','Bacillus subtilis subsp. subtilis str. 168','Mycobacterium tuberculosis H37Rv','Geobacter metallireducens GS-15','Clostridioides difficile 630']
taxonames = ['T. maritima MSB8','E. coli MG1655','P. putida KT2440','B. subtilis 168','M. tuberculosis H37Rv','G. metallireducens GS-15','C. difficile 630']
# taxonames = ['T. mari','E. coli','P. putida','B. sub','M. tuber','G. meta','C. diff']

data['newtotal'] = data['newtotal'].astype(float)
data['total_score'] = data['total_score'].astype(float)
data = data.sort_values(by=['modelname', 'evaltype', 'iteration'])
models = data['modelname'].unique()

#########################################################################
# # Create a figure with 4 subplots, one for each eval
# # List of evaluations to plot
# evals = ['total_score', 'newtotal','total_reactions','total_genes','total_metabolites', 'total_compartments','consistency', 'metabolic_coverage','annotation_rxn']

# # Iterate over each eval and corresponding subplot
# # 设置大图网格
# n_rows = (len(evals) + 2) // 3  # 每行最多放4个小图
# fig, axes = plt.subplots(n_rows, 3, figsize=(20, 5 * n_rows), squeeze=False)

# # 遍历每个模型，生成对应的小图
# for id, eval in enumerate(evals):
#     ax = axes[id // 3][id % 3]  # 选择子图
#     for idx, model in enumerate(models):
#         model_data = data[data['modelname'] == model]
#         model_data = model_data.sort_values(by='iteration')
#         # ax.plot(model_data['iteration'], model_data[eval], alpha=0.7)
#     ax.scatter(data['iteration'], data[eval],marker='.',alpha=0.7)
#     # data.scatter(column=[eval], by='iteration', ax=ax)
#     ax.set_ylim(data[eval].min() - 0.01, data[eval].max() + 0.01)
#     ax.set_title(eval)
#     ax.set_xlabel('Iteration')
#     ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
#     ax.set_ylabel(eval)
# # 移除多余的空白子图
# for j in range(len(models), n_rows * 3):
#     fig.delaxes(axes[j // 3][j % 3])
# plt.suptitle('')
# plt.tight_layout()
# plt.savefig(f'../data/result/figures/eval7_memote_t2up_boxplot.png', dpi=300)
# exit()


#########################################################################
import seaborn as sns
# print(data.columns)
print(data.head())

evals = ['total_reactions','total_genes','total_metabolites','consistency', 'metabolic_coverage']
evals = ['total_score','consistency', 'metabolic_coverage','total_reactions','total_genes','total_metabolites']

filtered_data = data[data['evaltype'] == 'differ']

# Melt the DataFrame into long format
melted_data = filtered_data.melt(
    id_vars=["iteration"],
    value_vars=evals,
    var_name="metric",
    value_name="value"
)

# Set up the FacetGrid
sns.set_theme(style="whitegrid")
g = sns.FacetGrid(melted_data, col="metric", col_wrap=3, height=4, sharey=False)

# Map a lineplot to the grid
g.map_dataframe(sns.lineplot, x="iteration", y="value", marker='o')

# Customize titles and labels
# subplot names are the metric names but without the "_" and Capitalized
colnames = [metric.replace("_", " ").capitalize() for metric in evals]
for ax, colname in zip(g.axes.flat, colnames):
    ax.set_title(colname)
    # ax.set_xlabel("Iteration")
    # ax.set_ylabel("Value")
g.set_axis_labels("Iteration", "Value")

# adjust y-axis limits by each metric
for ax, metric in zip(g.axes.flat, evals):
    metric_data = filtered_data[filtered_data['evaltype'] == 'differ']
    metric_data = metric_data[metric]
    if metric == 'consistency' or metric == 'total_score':
        ax.set_ylim(metric_data.min() - 0.1, metric_data.max() + 0.1)


# Adjust layout and show the plot

plt.tight_layout()
plt.savefig(f'../data/result/figures/eval7_memote_t2up_lineplot.png', dpi=300)

















n_rows = (len(evals) + 2) // 3  # 每行最多放4个小图
fig, axes = plt.subplots(n_rows, 3, figsize=(20, 5 * n_rows), squeeze=False)

# 遍历每个模型，生成对应的小图
for id, eval in enumerate(evals):
    ax = axes[id // 3][id % 3]  # 选择子图
    for idx, model in enumerate(models):
        model_data = data[data['modelname'] == model]
        model_data = model_data.sort_values(by='iteration')
        ax.plot(model_data['iteration'], model_data[eval], alpha=0.7)
    # ax.scatter(data['iteration'], data[eval],marker='.',alpha=0.7)
    # data.scatter(column=[eval], by='iteration', ax=ax)
    ax.set_ylim(data[eval].min() - 0.01, data[eval].max() + 0.01)
    ax.set_title(eval)
    ax.set_xlabel('Iteration')
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax.set_ylabel(eval)
# 移除多余的空白子图
for j in range(len(evals), n_rows * 3):
    fig.delaxes(axes[j // 3][j % 3])
plt.suptitle('')
plt.tight_layout()
plt.savefig(f'../data/result/figures/eval7_memote_t2up_creasing.png', dpi=300)