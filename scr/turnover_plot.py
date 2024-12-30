import matplotlib.pyplot as plt
import pandas as pd
f ='/ibex/user/niuk0a/funcarve/cobra/data/result/eval_turnover_7v4.csv'
f ='/ibex/user/niuk0a/funcarve/cobra/data/result/eval_turnover_7v4_I1.csv'

# modelnames=(iLJ478 iJN1463 iCN900 iAF1260 iAF987 iNJ661 iYO844)
# modelname = 'iJN1463'
# df = pd.read_csv(f)
# df = df.dropna().drop_duplicates()
# df.columns = ['modelname', 'types', 'reward', 'threshold', 'allpr_change','have_label','full_match','label_rate','full_match_rate','acc','f1','precision','recall']
# df = df[df['modelname']== modelname]
# df =df[df['types']=='allec']

# ## draw boxplot for f1, acc, recall, precision and the box contains the median, 25% and 75% quantile
# ## x-axis is f1, acc, recall, precision
# ## y-axis is the value of f1, acc, recall, precision


# metrics = ['f1', 'acc', 'recall', 'precision']
# data_to_plot = [df[metric] for metric in metrics]
# print(data_to_plot)
# plt.figure(figsize=(10, 6))
# plt.boxplot(data_to_plot, labels=metrics,)
# plt.xlabel('Metrics')
# plt.ylabel('Value')
# plt.title('Boxplot of F1, Accuracy, Recall, and Precision')
# plt.grid(True)
# plt.show()

# plt.savefig(f'/ibex/user/niuk0a/funcarve/cobra/data/result/figures/eval_turnover_7v4_{modelname}.png')
# exit()
#################
# # iLJ478	NC_000853.1	243274	Thermotoga maritima MSB8
# # iJN1463	NC_002947.4	160488	Pseudomonas putida KT2440
# # iCN900	AM180355.1	272563	Clostridioides difficile 630
# # iAF1260	NC_000913.3	83333	Escherichia coli str. K-12 substr. MG1655
# # iAF987	CP000148.1	269799	Geobacter metallireducens GS-15
# # iNJ661	NC_000962.3	83332	Mycobacterium tuberculosis H37Rv
# # iYO844	AL009126.3	224308	Bacillus subtilis subsp. subtilis str. 168
# modelnames = ['iLJ478', 'iJN1463', 'iCN900', 'iAF1260', 'iAF987', 'iNJ661', 'iYO844']
# taxonames=['Thermotoga maritima MSB8','Pseudomonas putida KT2440','Clostridioides difficile 630','Escherichia coli str. K-12 substr. MG1655','Geobacter metallireducens GS-15','Mycobacterium tuberculosis H37Rv','Bacillus subtilis subsp. subtilis str. 168']
# df = pd.read_csv(f)
# df = df.dropna().drop_duplicates()
# df.columns = ['modelname', 'types', 'reward', 'threshold', 'allpr_change','have_label','full_match','label_rate','full_match_rate','acc','f1','precision','recall']
# df = df[df['types']=='allec']


# metrics = ['f1', 'acc', 'recall', 'precision']
# fig, axes = plt.subplots(nrows=2, ncols=(len(modelnames) + 1) // 2, figsize=(14, 7), sharey=True)

# for ax, modelname, taxon in zip(axes.flat, modelnames, taxonames):
#     df_model = df[df['modelname'] == modelname]
#     print(modelname)
#     print(df_model)
#     data_to_plot = [df_model[metric] for metric in metrics]
#     ax.boxplot(data_to_plot, labels=metrics)
#     ax.set_title(f'{taxon}')
#     ax.set_xlabel('Metrics')
#     ax.grid(True)


# axes[0, 0].set_ylabel('Value')
# axes[1, 0].set_ylabel('Value')
# ## remove empty subplots
# for i in range(len(modelnames), len(axes.flat)):
#     fig.delaxes(axes.flat[i])
# fig.suptitle('Boxplot of F1, Accuracy, Recall, and Precision for Different Organisms')
# plt.tight_layout(rect=[0, 0.03, 1, 0.95])
# plt.show()

# fig.savefig('/ibex/user/niuk0a/funcarve/cobra/data/result/figures/eval_turnover_7v4_all_models.png')
# exit()

#########################################################################################
# only iter1
f = '/ibex/user/niuk0a/funcarve/cobra/data/result/eval_turnover_7v4_I1.csv'
df = pd.read_csv(f)
df = df.dropna().drop_duplicates()
df.columns = ['modelname', 'types', 'reward', 'threshold', 'allpr_change','have_label','full_match','label_rate','full_match_rate','acc','f1','precision','recall','below_t','below_t_all_match']
# modelnames = ['iLJ478', 'iJN1463', 'iCN900', 'iAF1260', 'iAF987', 'iNJ661', 'iYO844']
# taxonames=['Thermotoga maritima MSB8','Pseudomonas putida KT2440','Clostridioides difficile 630','Escherichia coli str. K-12 substr. MG1655','Geobacter metallireducens GS-15','Mycobacterium tuberculosis H37Rv','Bacillus subtilis subsp. subtilis str. 168']

modelnames = ['iLJ478','iAF1260','iJN1463','iYO844','iNJ661','iAF987','iCN900']
taxonames = ['Thermotoga maritima MSB8','Escherichia coli str. K-12 substr. MG1655','Pseudomonas putida KT2440','Bacillus subtilis subsp. subtilis str. 168','Mycobacterium tuberculosis H37Rv','Geobacter metallireducens GS-15','Clostridioides difficile 630']
taxonames = ['T. maritima MSB8','E. coli MG1655','P. putida KT2440','B. subtilis 168','M. tuberculosis H37Rv','G. metallireducens GS-15','C. difficile 630']
taxonames = ['T. mari','E. coli','P. putida','B. sub','M. tuber','G. meta','C. diff']
# shorten the taxonames
# taxonames
# the plot for all models and all metrics in one figure
metrics = ['f1', 'acc', 'recall', 'precision']
metricsname = ['F1', 'Accuracy', 'Recall', 'Precision']
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 5), sharey=True)
for ax, metric in zip(axes.flat, metrics):
    data_to_plot = [df[df['modelname'] == modelname][metric] for modelname in modelnames]
    ax.boxplot(data_to_plot, labels=taxonames)
    # ax.set_title(f'{metric}')
    ax.set_title(f'{metricsname[metrics.index(metric)]}')

    ax.set_xlabel('Genome')
    ax.grid(True)
    ## add a horizontal line at y=0.5
    ax.axhline(y=0.5, color='grey', linestyle='--')
axes[0, 0].set_ylabel('Value')
axes[1, 0].set_ylabel('Value')
fig.suptitle('Boxplot of F1, Accuracy, Recall, and Precision for Different Models', fontsize=16)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()
fig.savefig('/ibex/user/niuk0a/funcarve/cobra/data/result/figures/eval_turnover_7v4_I1_all_models.png')

