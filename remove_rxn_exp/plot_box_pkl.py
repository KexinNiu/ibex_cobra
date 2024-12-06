from builtins import zip
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pickle


f ='/ibex/user/niuk0a/funcarve/cobra/remove_rxn_exp/01result.csv'
f ='/ibex/user/niuk0a/funcarve/cobra/remove_rxn_exp/01result_cleandf.csv'
# f='/ibex/user/niuk0a/funcarve/cobra/remove_rxn_exp/01t.csv'
orf='/ibex/user/niuk0a/funcarve/cobra/remove_rxn_exp/01result_ori.csv'


## read the data
df = pd.read_csv(f,sep='\t',skiprows=1,header=None)
df.columns = ['prec','seednum','rm_rxn','rm_rxn_id','method','ori_allgfreactions','wgf_allgfreactions','ori_tp','ori_fn','ori_fp','ori_f1','wgf_tp','wgf_fn','wgf_fp','wgf_f1']


# remove duplicate rows
df = df.drop_duplicates()
df['prec'] = df['prec'].astype(float)
df['seednum'] = df['seednum'].astype(int)
df['rm_rxn_id'] = df['rm_rxn_id'].apply(lambda x: x[1:-1].split(','))
df['wgf_allgfreactions'] = df['wgf_allgfreactions'].apply(lambda x: int(x[1:-1]))
df['wgf_tp'] = df['wgf_tp'].apply(lambda x: int(x[1:-1]))
df['wgf_fn'] = df['wgf_fn'].apply(lambda x:int(x[1:-1]))
df['wgf_fp'] = df['wgf_fp'].apply(lambda x: int(x[1:-1]))
df['wgf_f1'] = df['wgf_f1'].apply(lambda x: float(x[1:-1]))
df['method'] = df['method'].apply(lambda x: int(x[1:-1]))
df['ori_allgfreactions'] = [0] * len(df)
df['ori_tp'] = [0] * len(df)
df['ori_fn'] = [0] * len(df)
df['ori_fp'] = [0] * len(df)
df['ori_f1'] = [0.0] * len(df)

print(df.head())
# prec  seednum                                          rm_rxn_id  method ori_allgfreactions  ...  ori_f1 wgf_tp wgf_fn wgf_fp    wgf_f1
# 0  0.05     4444  ['MLDEP2pp',  'SADT2',  'GLUSx',  'CRO4tex',  ...       1                 []  ...      []      4     15    736  0.010540
# 1  0.05     4444  ['MLDEP2pp',  'SADT2',  'GLUSx',  'CRO4tex',  ...       2                 []  ...      []      1     18    144  0.012195
# 2  0.05     6666  ['CITCIb',  'AGMH',  'Cobalt2abcppI',  'MCTP1B...       1                 []  ...      []      3     16    743  0.007843
# 3  0.05     4444  ['MLDEP2pp',  'SADT2',  'GLUSx',  'CRO4tex',  ...       3                 []  ...      []      1     18     13  0.060606

# [4 rows x 14 columns]

## read the data for original
df_ori = pd.read_csv(orf,sep='\t',header=None)
df_ori.columns = ['prec','seednum','rm_rxn','rm_rxn_id','method','ori_allgfreactions','wgf_allgfreactions','ori_tp','ori_fn','ori_fp','ori_f1','wgf_tp','wgf_fn','wgf_fp','wgf_f1']
df_ori['prec'] = df_ori['prec'].astype(float)
df_ori['seednum'] = df_ori['seednum'].astype(int)
df_ori['rm_rxn_id'] = df_ori['rm_rxn_id'].apply(lambda x: x[1:-1].split(','))
df_ori['ori_allgfreactions'] = df_ori['ori_allgfreactions'].apply(lambda x: int(x[1:-1]))
df_ori['ori_tp'] = df_ori['ori_tp'].apply(lambda x: int(x[1:-1]))
df_ori['ori_fn'] = df_ori['ori_fn'].apply(lambda x:int(x[1:-1]))
df_ori['ori_fp'] = df_ori['ori_fp'].apply(lambda x: int(x[1:-1]))
df_ori['ori_f1'] = df_ori['ori_f1'].apply(lambda x: float(x[1:-1]))
print(df_ori.head())

#    prec  seednum  rm_rxn                                          rm_rxn_id method  ...    ori_f1 wgf_tp  wgf_fn  wgf_fp  wgf_f1
# 0  0.50     6666     197  ['CITCIb',  'AGMH',  'Cobalt2abcppI',  'MCTP1B...     []  ...  0.009346     []      []      []      []
# 1  0.25     6666      98  ['CITCIb',  'AGMH',  'Cobalt2abcppI',  'MCTP1B...     []  ...  0.000000     []      []      []      []
# 2  0.40     6666     157  ['CITCIb',  'AGMH',  'Cobalt2abcppI',  'MCTP1B...     []  ...  0.000000     []      []      []      []
# 3  0.50     5555     197  ['GCCa',  'ACtex',  '4HBAtex',  'COBALT2tex', ...     []  ...  0.039344     []      []      []      []
# 4  0.45     6666     177  ['CITCIb',  'AGMH',  'Cobalt2abcppI',  'MCTP1B...     []  ...  0.010363     []      []      []      []

# [5 rows x 15 columns]


## merge the data by seednum 
print('-------------------')

## we aim to draw the boxplot for the data which consist of the following :1.ori_f1 2.wgf_f1  x-axis: prec, y-axis: f1 and the boxplot need to consider different seednum

# Merge the dataframes to update 'df' with values from 'df_ori'
# Ensure consistent data types in both dataframes
columns_to_update = ['ori_allgfreactions', 'ori_tp', 'ori_fn', 'ori_fp', 'ori_f1']

for col in columns_to_update:
    # Convert both dataframes' columns to a common dtype
    df[col] = pd.to_numeric(df[col], errors='coerce')
    df_ori[col] = pd.to_numeric(df_ori[col], errors='coerce')

# Merge the dataframes on 'prec' and 'seednum', filling missing ori_* values in df from df_ori
df = df.set_index(['prec', 'seednum'])
df_ori = df_ori.set_index(['prec', 'seednum'])
df.update(df_ori[columns_to_update])

# Reset the index after updating
df = df.reset_index()

# Verify the updated dataframe
print(df.head())
# Verify the updated dataframe by checking the first few rows
print(df[['prec', 'seednum','method', 'ori_f1', 'wgf_f1','rm_rxn']].head())
print(df[df['prec'] == 0.5].head())
# Prepare data for plotting
#########################################F1 SCORE METRICS#############################################333
# plot_data = df[['prec', 'seednum', 'method', 'ori_f1', 'wgf_f1','wgf_allgfreactions','ori_allgfreactions']]
plot_data = df[['prec', 'seednum', 'method', 'ori_f1', 'wgf_f1']]

# Melt the dataframe for easier plotting
plot_data = plot_data.melt(id_vars=['prec', 'seednum', 'method'], 
                           var_name='Metric', 
                           value_name='F1 Score')

# Create the boxplot, adding 'method' as a hue parameter
plt.figure(figsize=(14, 7))
sns.boxplot(data=plot_data, x='prec', y='F1 Score', hue='Metric', palette="Set2")

# Add facet labels for different methods
g = sns.catplot(data=plot_data, 
                x='prec', 
                y='F1 Score', 
                hue='Metric', 
                col='method', 
                kind='box', 
                col_wrap=3, 
                height=4, 
                aspect=1.2, 
                palette='Set2')

# Configure the plots
g.set_titles("Method {col_name}")
g.set_axis_labels("Precision (prec)", "F1 Score")
g.tight_layout()

# Show the plot
plt.show()
plt.savefig('/ibex/user/niuk0a/funcarve/cobra/remove_rxn_exp/plot_test1_cleandf.png')
exit()


#### draw another plot for tp 
plot_data = df[['prec', 'seednum', 'method', 'ori_tp', 'wgf_tp','rm_rxn']]

# Prepare data for boxplots
boxplot_data = df[['prec', 'seednum', 'method', 'ori_tp', 'wgf_tp']].melt(
    id_vars=['prec', 'seednum', 'method'], 
    value_name='TP', 
    var_name='Metric'
)

# Prepare data for line plot (rm_rxn)
lineplot_data = df[['prec', 'method', 'rm_rxn']]

# Create a FacetGrid for subplots by method
g = sns.FacetGrid(boxplot_data, col='method', col_wrap=3, height=4, aspect=1.2)

# Add boxplots for `ori_tp` and `wgf_tp` to each subplot
g.map_dataframe(
    sns.boxplot, 
    x='prec', 
    y='TP', 
    hue='Metric', 
    data=boxplot_data,
    palette='Set2'
)

# Overlay line plots for `rm_rxn` on the same subplots
for ax, (method, data) in zip(g.axes.flat, lineplot_data.groupby('method')):
    # Convert the categorical x-axis to numerical for alignment
    prec_categories = sorted(boxplot_data['prec'].unique())
    prec_mapping = {cat: idx for idx, cat in enumerate(prec_categories)}
    
    # Map `prec` to its categorical positions
    data['prec_numeric'] = data['prec'].map(prec_mapping)
    
    sns.lineplot(
        data=data, 
        x='prec_numeric', 
        y='rm_rxn', 
        marker='o', 
        ax=ax, 
        color='red', 
        label='rm_rxn'
    )
    
    # Update x-tick labels to match the categorical `prec`
    ax.set_xticks(range(len(prec_categories)))
    ax.set_xticklabels(prec_categories)
# Adjust legends and titles
for ax in g.axes.flat:
    handles, labels = ax.get_legend_handles_labels()
    if "rm_rxn" in labels:  # Combine boxplot and lineplot legends
        ax.legend(loc='upper right')
## custom the y-axis range
    ax.set_ylim(0, 30)
g.set_titles("Method: {col_name}")
g.set_axis_labels("Precision (prec)", "TP")
g.tight_layout()

# Save and show the plot
plt.savefig('/ibex/user/niuk0a/funcarve/cobra/remove_rxn_exp/plot_tp1_cleandf.png')
plt.show()




exit()


###############################################################################3
# Prepare data for plotting

# bf = open(f,'rb')
# dd = pickle.load(bf)

# df = pd.DataFrame(dd)

# # Create a DataFrame from the provided data

# df = pd.DataFrame(data)

# # resultdd['prec'].append(result['prec'][i])
# # resultdd['seednum'].append(result['seednum'][i])
# # resultdd['rm_rxn_id'].append(result['rm_rxn_id'][i])
# # resultdd['ori_allgfreactions'].append(result['ori_allgfreactions'][i])
# # resultdd['wgf_allgfreactions'].append(result['wgf_allgfreactions'][i])
# # resultdd['ori_tp'].append(result['ori_tp'][i])
# # resultdd['ori_fn'].append(result['ori_fn'][i])
# # resultdd['ori_fp'].append(result['ori_fp'][i])
# # resultdd['wgf_tp'].append(result['wgf_tp'][i])
# # resultdd['wgf_fn'].append(result['wgf_fn'][i])
# # resultdd['wgf_fp'].append(result['wgf_fp'][i])
# # data = data_unre


# # Prepare the data for boxplot
# addW_reactions_flat = pd.DataFrame(df['addW_reactions'].tolist(), index=df['percentage']).melt(var_name="addW_index", value_name="addW_reactions", ignore_index=False).reset_index()
# addW_reactions_inRM_num_flat = pd.DataFrame(df['addW_reactions_inRM_num'].tolist(), index=df['percentage']).melt(var_name="inRM_index", value_name="addW_reactions_inRM_num", ignore_index=False).reset_index()
# addgf_reactions_inRM_num_flat = pd.DataFrame(df['addgf_reactions_inRM_num'].tolist(), index=df['percentage']).melt(var_name="inRM_index", value_name="addgf_reactions_inRM_num", ignore_index=False).reset_index()

# wgf_fp = pd.DataFrame(resultdd['wgf_fp']).melt(var_name="inRM_index", value_name="addgf_reactions_inRM_num", ignore_index=False).reset_index()

# # Set up the plot
# fig, ax = plt.subplots(figsize=(12, 8))

# # Convert 'percentage' to categorical for proper spacing
# addW_reactions_flat['percentage'] = addW_reactions_flat['percentage'].astype(str)
# addW_reactions_inRM_num_flat['percentage'] = addW_reactions_inRM_num_flat['percentage'].astype(str)
# df['percentage'] = df['percentage'].astype(str)
# ################# for box plot #################二选一
# # Boxplot for addW_reactions
# box1 = sns.boxplot(x='percentage', y='addW_reactions', data=addW_reactions_flat, ax=ax, color='Yellow', boxprops=dict(alpha=0.5), label='W_add_reactions')

# # Boxplot for addW_reactions_inRM_num
# box2 =sns.boxplot(x='percentage', y='addW_reactions_inRM_num', data=addW_reactions_inRM_num_flat, ax=ax, color='Green', boxprops=dict(alpha=0.5),label='W_correct_reactions')
# box3 = sns.boxplot(x='percentage', y='addgf_reactions_inRM_num', data=addgf_reactions_inRM_num_flat, ax=ax, color='grey', boxprops=dict(alpha=0.5),label='gf_correct_reactions')

# # Plot for rm_reactions using sns.lineplot to connect the points with the correct spacing
# linp = sns.lineplot(x='percentage', y='rm_reactions', data=df, marker='o', color='darkgrey',label='Remove reactions', ax=ax, linewidth=1, markersize=3)

# # Customize the plot
# ax.set_title('Boxplot of Weighted/Non-weighted gapfilling(All unreviewed genes)')
# ax.set_xlabel('Percentage')
# ax.set_ylabel('Reaction Counts')

# # Add a legend
# ax.legend(title='Legend', labels=['Remove reactions', 'W_add_reactions', 'W_correct_reactions', 'gf_correct_reactions'])

# plt.legend()
# plt.show()
# # plt.savefig('boxplot.png')
# plt.savefig('boxplot_allun.png')
# ################# for box plot #################

# ################# for compare #################二选一
# # addW_reactions_inRM_num_flat = pd.DataFrame(df['addW_reactions_inRM_num'].tolist(), index=df['percentage']).melt(var_name="inRM_index", value_name="addW_reactions_inRM_num", ignore_index=False).reset_index()
# # addW_reactions_inRM_num_unre_flat = pd.DataFrame(df['addW_reactions_inRM_num_unre'].tolist(), index=df['percentage']).melt(var_name="inRM_index", value_name="addW_reactions_inRM_num_unre", ignore_index=False).reset_index()
# # boxcorr =sns.boxplot(x='percentage', y='addW_reactions_inRM_num', data=addW_reactions_inRM_num_flat, ax=ax, color='Yellow', boxprops=dict(alpha=0.5),label='W_correct_reactions')
# # boxcorrun =sns.boxplot(x='percentage', y='addW_reactions_inRM_num_unre', data=addW_reactions_inRM_num_unre_flat, ax=ax, color='Green', boxprops=dict(alpha=0.5),label='W_correct_reactions_unreviewed')

# # ax.set_title('Boxplot of Weighted gapfilling Correct reactions Have/Only-have Unreviewed genes)')
# # ax.set_xlabel('Percentage')
# # ax.set_ylabel('Reaction Counts')
# # ax.legend(title='Legend', labels=['Remove reactions', 'W_correct_reactions', 'W_correct_reactions_unreviewed'])
# # plt.legend()

# # plt.savefig('boxplot_compare.png')
################# for compare #################