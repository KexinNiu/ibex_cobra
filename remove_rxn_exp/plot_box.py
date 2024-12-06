import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Create a DataFrame from the provided data
data = {
    'percentage': [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5],
    'addW_reactions': [[37, 37, 34, 38, 22], [38, 39, 39, 47, 24], [44, 57, 51, 48, 27],
                       [66, 67, 76, 53, 46], [74, 71, 108, 69, 40], [71, 75, 111, 76, 47],
                       [73, 75, 118, 82, 60], [83, 81, 125, 94, 69], [96, 84, 123, 96, 76],
                       [98, 85, 124, 103, 78]],
    'rm_reactions': [21,43,65,87,109,131,153,175,197,219],
    # 'addgf_reactions': [[37, 37, 34, 38, 22], [38, 39, 39, 47, 24], [44, 57, 51, 48, 27],
    #                    [66, 67, 76, 53, 46], [74, 71, 108, 69, 40], [71, 75, 111, 76, 47],
    #                    [73, 75, 118, 82, 60], [83, 81, 125, 94, 69], [96, 84, 123, 96, 76],
    #                    [98, 85, 124, 103, 78]],
    
    'addgf_reactions_inRM_num': [[0,0,0,0,0], [0,0,0,0,0], [0,0,0,0,0],[0,0,0,0,0],
                                 [0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],
                                 [0,0,0,0,0],[0,0,0,0,0]],
    'addW_reactions_inRM_num': [[3,3, 2, 1, 1], [5,3, 3, 3, 2], [6,2, 7, 4, 2], [5,8, 9, 6, 4],
                                [5,9, 11, 4, 9], [5,11, 14, 6, 8], [8,13, 18, 9, 10], [13,13, 17, 16, 10],
                                [15,14, 18, 18, 12], [19,15, 19, 21, 17]],
                                'addW_reactions_inRM_num_unre': [
[2,1, 1,2,3], 
[5,3, 1, 4, 7], 
[7,3, 5, 9, 10], 
[9,8, 6, 11, 8],
[7,10, 8, 7, 11], 
[12,9, 10, 11, 13], 
[6,9, 13, 13, 9], 
[14,10, 12, 17, 11],
[12,11, 11, 17, 13], 
[17,13, 14, 15, 15]]
    
}
data_unre ={
    'percentage': [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5],
    'addW_reactions': [
[27,29,28,27,30],
[35,35,31,31,38],
[39,41,43,39,47],
[47,48,50,44,52],
[51,55,48,57,59],
[54,54,52,66,63],
[53,72,61,72,67],
[62,72,78,93,69],
[64,84,81,89,75],
[73,93,87,89,80]
],
    'rm_reactions': [19,39,59,78,98,118,137,157,177,197],
    'addgf_reactions_inRM_num': [[0,0,0,0,0], [0,0,0,0,0], [0,0,0,0,0],[0,0,0,0,0],
                                 [0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],
                                 [0,0,0,0,0],[0,0,0,0,0]],
    'addW_reactions_inRM_num': [
[2,1, 1,2,3], 
[5,3, 1, 4, 7], 
[7,3, 5, 9, 10], 
[9,8, 6, 11, 8],
[7,10, 8, 7, 11], 
[12,9, 10, 11, 13], 
[6,9, 13, 13, 9], 
[14,10, 12, 17, 11],
[12,11, 11, 17, 13], 
[17,13, 14, 15, 15]]

}

# data = data_unre
df = pd.DataFrame(data)

# Prepare the data for boxplot
addW_reactions_flat = pd.DataFrame(df['addW_reactions'].tolist(), index=df['percentage']).melt(var_name="addW_index", value_name="addW_reactions", ignore_index=False).reset_index()
addW_reactions_inRM_num_flat = pd.DataFrame(df['addW_reactions_inRM_num'].tolist(), index=df['percentage']).melt(var_name="inRM_index", value_name="addW_reactions_inRM_num", ignore_index=False).reset_index()
addgf_reactions_inRM_num_flat = pd.DataFrame(df['addgf_reactions_inRM_num'].tolist(), index=df['percentage']).melt(var_name="inRM_index", value_name="addgf_reactions_inRM_num", ignore_index=False).reset_index()

# Set up the plot
fig, ax = plt.subplots(figsize=(12, 8))

# Convert 'percentage' to categorical for proper spacing
addW_reactions_flat['percentage'] = addW_reactions_flat['percentage'].astype(str)
addW_reactions_inRM_num_flat['percentage'] = addW_reactions_inRM_num_flat['percentage'].astype(str)
df['percentage'] = df['percentage'].astype(str)
################# for box plot #################二选一
# Boxplot for addW_reactions
box1 = sns.boxplot(x='percentage', y='addW_reactions', data=addW_reactions_flat, ax=ax, color='Yellow', boxprops=dict(alpha=0.5), label='W_add_reactions')

# Boxplot for addW_reactions_inRM_num
box2 =sns.boxplot(x='percentage', y='addW_reactions_inRM_num', data=addW_reactions_inRM_num_flat, ax=ax, color='Green', boxprops=dict(alpha=0.5),label='W_correct_reactions')
box3 = sns.boxplot(x='percentage', y='addgf_reactions_inRM_num', data=addgf_reactions_inRM_num_flat, ax=ax, color='grey', boxprops=dict(alpha=0.5),label='gf_correct_reactions')

# Plot for rm_reactions using sns.lineplot to connect the points with the correct spacing
linp = sns.lineplot(x='percentage', y='rm_reactions', data=df, marker='o', color='darkgrey',label='Remove reactions', ax=ax, linewidth=1, markersize=3)

# Customize the plot
ax.set_title('Boxplot of Weighted/Non-weighted gapfilling(All unreviewed genes)')
ax.set_xlabel('Percentage')
ax.set_ylabel('Reaction Counts')

# Add a legend
ax.legend(title='Legend', labels=['Remove reactions', 'W_add_reactions', 'W_correct_reactions', 'gf_correct_reactions'])

plt.legend()
plt.show()
plt.savefig('boxplot.png')
plt.savefig('boxplot_allun.png')
################# for box plot #################

################# for compare #################二选一
# addW_reactions_inRM_num_flat = pd.DataFrame(df['addW_reactions_inRM_num'].tolist(), index=df['percentage']).melt(var_name="inRM_index", value_name="addW_reactions_inRM_num", ignore_index=False).reset_index()
# addW_reactions_inRM_num_unre_flat = pd.DataFrame(df['addW_reactions_inRM_num_unre'].tolist(), index=df['percentage']).melt(var_name="inRM_index", value_name="addW_reactions_inRM_num_unre", ignore_index=False).reset_index()
# boxcorr =sns.boxplot(x='percentage', y='addW_reactions_inRM_num', data=addW_reactions_inRM_num_flat, ax=ax, color='Yellow', boxprops=dict(alpha=0.5),label='W_correct_reactions')
# boxcorrun =sns.boxplot(x='percentage', y='addW_reactions_inRM_num_unre', data=addW_reactions_inRM_num_unre_flat, ax=ax, color='Green', boxprops=dict(alpha=0.5),label='W_correct_reactions_unreviewed')

# ax.set_title('Boxplot of Weighted gapfilling Correct reactions Have/Only-have Unreviewed genes)')
# ax.set_xlabel('Percentage')
# ax.set_ylabel('Reaction Counts')
# ax.legend(title='Legend', labels=['Remove reactions', 'W_correct_reactions', 'W_correct_reactions_unreviewed'])
# plt.legend()

# plt.savefig('boxplot_compare.png')
################# for compare #################