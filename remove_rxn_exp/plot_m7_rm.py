import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
pref='/ibex/user/niuk0a/funcarve/cobra/remove_rxn_exp/mmodel7_iAF987_ori.csv'
modelnames = [
    'iLJ478',
    'iJN1463',#1
    'iCN900',
    'iAF1260',#3
    'iAF987',
    'iNJ661',
    'iYO844'#6
]
# modelnames = modelnames[1,3,6]
modelnames = ['iJN1463','iAF1260','iYO844']
modelnames = ['iAF987']
#iLJ478	Thermotoga maritima MSB8
# iJN1463	Pseudomonas putida KT2440
# iCN900	Clostridioides difficile 630
# iAF1260	Escherichia coli str. K-12 substr. MG1655
# iAF987	Geobacter metallireducens GS-15
# iNJ661	Mycobacterium tuberculosis H37Rv
# iYO844	Bacillus subtilis subsp. subtilis str. 168
taxonames=[
    'Thermotoga maritima MSB8',
    'Pseudomonas putida KT2440',
    'Clostridioides difficile 630',
    'Escherichia coli str. K-12 substr. MG1655',
    'Geobacter metallireducens GS-15',
    'Mycobacterium tuberculosis H37Rv',
    'Bacillus subtilis subsp. subtilis str. 168'
]
# taxonames = [taxonames[1],taxonames[3],taxonames[4],taxonames[6]]
taxonames = [taxonames[4]]
typenames = ['ori','wgf']
evals = ['f1','acc','coverage','addrxncount']
# evals = ['f1','coverage']
# evals = ['recall','precision']
# evals = ['addrxncount']

for model in modelnames:
    for evalnow in evals:
        print('current model:', model,'evalus:',evalnow)
        orif = f'/ibex/user/niuk0a/funcarve/cobra/remove_rxn_exp/mmodel7_{model}_{typenames[0]}_r.csv'
        oridf = pd.read_csv(orif,sep='\t')
        # print(oridf.head())
        oridf.columns =['modelname','prec','seednum','rm_rxn_id_nub','rm_rxn_id','method','ori_addrxncount',
                        'ori_tp','ori_fn','ori_fp','ori_f1']
        wgf = f'/ibex/user/niuk0a/funcarve/cobra/remove_rxn_exp/mmodel7_{model}_{typenames[1]}_r.csv'
        wgdf = pd.read_csv(wgf,sep='\t')
        # print(wgdf.head())
        wgdf.columns =['modelname','prec','seednum','rm_rxn_id_nub','rm_rxn_id','method','wgf_addrxncount',
                        'wgf_tp','wgf_fn','wgf_fp','wgf_f1']
        
        # evalmetric='coverage'
        evalmetric = evalnow
        ## format the data to float and int
        oridf['ori_f1'] = oridf['ori_f1'].apply(lambda x: float(x[1:-1]))
        oridf['ori_tp'] = oridf['ori_tp'].apply(lambda x: int(x[1:-1]))
        oridf['ori_fn'] = oridf['ori_fn'].apply(lambda x: int(x[1:-1]))
        oridf['ori_fp'] = oridf['ori_fp'].apply(lambda x: int(x[1:-1]))
        oridf['rm_rxn_id_nub'] = oridf['rm_rxn_id_nub'].astype(int)
        oridf['ori_addrxncount'] = oridf['ori_addrxncount'].apply(lambda x: int(x[1:-1]))
        oridf['prec'] = oridf['prec'].astype(float)
        wgdf['wgf_f1'] = wgdf['wgf_f1'].apply(lambda x: float(x[1:-1]))
        wgdf['wgf_tp'] = wgdf['wgf_tp'].apply(lambda x: int(x[1:-1]))
        wgdf['wgf_fn'] = wgdf['wgf_fn'].apply(lambda x: int(x[1:-1]))
        wgdf['wgf_fp'] = wgdf['wgf_fp'].apply(lambda x: int(x[1:-1]))
        wgdf['prec'] = wgdf['prec'].astype(float)
        wgdf['method'] = wgdf['method'].apply(lambda x: int(x[1:-1]))
        wgdf['rm_rxn_id_nub'] = wgdf['rm_rxn_id_nub'].astype(int)
        wgdf['wgf_addrxncount'] = wgdf['wgf_addrxncount'].apply(lambda x: int(x[1:-1]))
        #merge oridf and wgdf on 'modelname','prec' and 'seednum'
        oridf = oridf[['modelname','prec','seednum','rm_rxn_id_nub','ori_tp','ori_fn','ori_fp','ori_f1','ori_addrxncount']]
        wgdf = wgdf[['modelname','prec','seednum','method','rm_rxn_id_nub','wgf_tp','wgf_fn','wgf_fp','wgf_f1','wgf_addrxncount']]
        # print(oridf.tail())
        # print(wgdf.tail())
        # wgdf = wgdf[wgdf['method']==1]
        df = pd.merge(oridf,wgdf,on=['modelname','prec','seednum','rm_rxn_id_nub'])
        # print(df.tail())
        
        if evalmetric =='f1':
            vn = 'F1 Score'
        elif evalmetric =='acc':
            # print('acc:')
            vn = 'Accuracy'
            df['ori_acc'] = df['ori_tp'] / df['ori_addrxncount']
            df['wgf_acc'] = df['wgf_tp']/ df['wgf_addrxncount']
            # print(df.head())
        elif evalmetric =='recall':
            vn = 'Recall'
            df['ori_recall'] = df['ori_tp'] / (df['ori_tp'] + df['ori_fn'])
            df['wgf_recall'] = df['wgf_tp'] / (df['wgf_tp'] + df['wgf_fn'])
        elif evalmetric =='precision':
            vn = 'Precision'
            df['ori_precision'] = df['ori_tp'] / (df['ori_tp'] + df['ori_fp'])
            df['wgf_precision'] = df['wgf_tp'] / (df['wgf_tp'] + df['wgf_fp'])
        elif evalmetric =='coverage':
            vn = 'Coverage'
            df['ori_coverage'] = df['ori_tp'] / df['rm_rxn_id_nub']
            df['wgf_coverage'] = df['wgf_tp'] / df['rm_rxn_id_nub']
        elif evalmetric =='addrxncount':
            vn = 'Added Reactions Count'
            

        # print(df.head())
        plot_data = df[['prec','seednum','method',f'ori_{evalmetric}',f'wgf_{evalmetric}']]
        # order the data
        plot_data = plot_data.sort_values(by=['prec','seednum','method'])

        plot_data = plot_data.melt(
            id_vars=['prec','seednum','method'],
            var_name='Metric',
            value_name=vn
        )
        # print(plot_data)
        # plt.figure(figsize=(14, 7))
        plt.figure(figsize=(5,10))
        sns.boxplot(data=plot_data, x='prec', y=vn, hue='Metric', palette="Set2")
        # Add facet labels for different methods
        g = sns.catplot(data=plot_data, 
                        x='prec', 
                        y=vn, 
                        hue='Metric', 
                        col='method', 
                        kind='box', 
                        col_wrap=3, 
                        # col_wrap=1, 
                        height=5, 
                        # aspect=1.2, 
                        aspect=1, 
                        palette='Set2')
        # sns.move_legend(g, "upper left", bbox_to_anchor=(.80, .25), frameon=True)
        sns.move_legend(g, "upper right", bbox_to_anchor=(.85, .95),frameon=True)
        # sns.move_legend
        ## set text size of legend items
        # g._legend.set_title('Metric',prop={'size':18})
        g._legend.set_title('Metric')

        g._legend.texts[0].set_text('Original Gapfilling')
        # g._legend.texts[0].set_fontsize(15)
        g._legend.texts[1].set_text('Weighted Gapfilling')
        ## move legend to no overlap with the plot
        g._legend.set_bbox_to_anchor((1, 0.5))
        # g._legend.texts[1].set_fontsize(15)
        
        g.set_titles('Settings: {col_name}',size=15)
        
         # g.set_titles('B. subtilis 168',size=15)
        # g.set_titles(taxonames[modelnames.index(model)],size=15)
        g.set_xlabels('Precetange of Reactions Removed',size=12)
        g.set_ylabels(vn,size=12)
        # g.set_xticklabels(size=12)
        # g.set_yticklabels(size=12)
        
        # g._legend.remove()
        
        # set one legend position at lower right

        # plt.legend(loc='lower right')
        # add title for the whole plot
        plt.suptitle(f'{taxonames[modelnames.index(model)]}', size=16)
        
        plt.tight_layout()
        # plt.savefig(f'/ibex/user/niuk0a/funcarve/cobra/remove_rxn_exp/mm{model}_ori_wgf_{evalmetric}.png')
        plt.savefig(f'/ibex/user/niuk0a/funcarve/cobra/remove_rxn_exp/EXP2_{model}_{evalmetric}.png')

    # break
    # break
