import pickle
import os
import pandas as pd
import cobra
import matplotlib.pyplot as plt
import seaborn as sns

prefs =[
    '/ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1',
    '/ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2',
    '/ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2',
    '/ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3_t4',
    '/ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_t4',
    '/ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_t5',
    '/ibex/user/niuk0a/CLEAN/app/results/inputs/AL009126.3_t4'
]
prefs = [i.replace('/ibex/user/niuk0a/CLEAN/app/results/inputs/','/ibex/user/niuk0a/funcarve/cobra/data/result/sbmls/') for i in prefs]
modelnames = [
    'iLJ478',
    'iJN1463',
    'iCN900',
    'iAF1260',
    'iAF987',
    'iNJ661',
    'iYO844'
]
latestmodelnames = [
    'iLJ478',
    'iJN1463',
    'iCN900',
    'iML1515',
    'iAF987',
    'iEK1008',
    'iYO844'
]
# thrlist=['p001','p005','p01','p03','p05','p07','1','2','3','5']
thrlist=['p02','p05','1','3']

def check_modelbigg(evaltype='allec'):

    with open('/ibex/user/niuk0a/funcarve/cobra/uniprot/bigg2mn.pkl', 'rb') as f:
        bigg2mn = pickle.load(f)
    with open('/ibex/user/niuk0a/funcarve/cobra/uniprot/mn2seed.pkl', 'rb') as f:
        mn2seed = pickle.load(f)
   
    # outf ='/ibex/user/niuk0a/funcarve/cobra/data/result/emodel_bigg_overlap_v4.csv'
    # if os.path.exists(outf):
    #     overlap_df = pd.read_csv(outf)
    # else:
    #     overlap_df = pd.DataFrame()
        
    overlap={
        'modelname':[],
        'thr':[],
        'threshod':[],
        'evaltype':[],
        'iter':[],
        'model_subtype':[],
        'both':[],
        'predmodel':[],
        'biggmodel':[],
        'ratio_predmodel':[],
        'ratio_biggmodel':[],
        'tp':[],
        'fp':[],
        'fn':[],
        'precision':[],
        'recall':[],
        'f1':[]
    }
    
    for i in range(0,len(modelnames)):
        modelname = modelnames[i]
        currn = latestmodelnames[i]
        
        biggmodelf = '/ibex/user/niuk0a/funcarve/cobra/uniprot/'+modelname+'.xml'
        
        biggmodel = cobra.io.read_sbml_model(biggmodelf)
        bigg_rxns = set(rxn.id for rxn in biggmodel.reactions)
        # bigg_mets = set(met.id for met in biggmodel.metabolites)
        # bigg_genes = set(gene.id for gene in biggmodel.genes)
        if currn != modelname:
            biggmodelfc = '/ibex/user/niuk0a/funcarve/cobra/uniprot/'+currn+'.xml'
            biggmodelc = cobra.io.read_sbml_model(biggmodelfc)
            bigg_rxnsc = set(rxn.id for rxn in biggmodelc.reactions)
            bigg_rxns = bigg_rxns | bigg_rxnsc

        mnrxns = set(bigg2mn[rxn] for rxn in bigg_rxns if rxn in bigg2mn.keys())
        seedrxns = set()
        c=0
        for mn in mnrxns:
            if mn not in mn2seed.keys():
                # print(mn)c
                c+=1
                pass
            else:
                for seed in mn2seed[mn]:
                    seedrxns.add(seed)
                
        print(c/len(mnrxns),'no mapping seed')
        print(len(bigg_rxns),len(mnrxns),len(seedrxns))
        print(len(seedrxns),list(seedrxns)[:10])
        ######
        reconmodel = cobra.io.read_sbml_model('/ibex/user/niuk0a/funcarve/cobra/uniprot/recon/'+modelname+'recon.sbml')
        recon_rxns = set(rxn.id.replace('_c','') for rxn in reconmodel.reactions)
        recon_mets = set(met.id for met in reconmodel.metabolites)
        recon_genes = set(gene.id for gene in reconmodel.genes)
        print('recon:',len(recon_rxns),len(recon_mets),len(recon_genes))
        print('recon rxn:',list(recon_rxns)[:10])
        ######
        intersection = len(recon_rxns & seedrxns)
        print('intersection:',intersection)  
        overlap['modelname'].append(modelname)
        overlap['thr'].append(0)
        overlap['threshod'].append('0')
        overlap['evaltype'].append('recon')
        overlap['iter'].append(0)
        overlap['model_subtype'].append('rxns')
        overlap['both'].append(intersection)
        overlap['predmodel'].append(len(recon_rxns))
        overlap['biggmodel'].append(len(seedrxns))
        overlap['ratio_predmodel'].append(intersection/len(recon_rxns))
        overlap['ratio_biggmodel'].append(intersection/len(seedrxns))
        overlap['tp'].append(intersection)
        overlap['fp'].append(len(recon_rxns)-intersection)
        overlap['fn'].append(len(seedrxns)-intersection)
        overlap['precision'].append(intersection/len(recon_rxns))
        overlap['recall'].append(intersection/len(seedrxns))
        overlap['f1'].append(2*intersection/(len(recon_rxns)+len(seedrxns)))
        # continue

  
        pref = prefs[i]
        # evaltype = 'allec'
        for thr in thrlist:
            for threshodl in ['5','6']:
                # 'NC_000853.1_t1_maxsep_dfiLJ478_allec_p01iter_1.sbml'
                print('>>',modelname,thr,flush=True)
                # for j in range(0,10):
                # /ibex/user/niuk0a/funcarve/cobra/data/result/sbmls/NC_002947.4_t2iJN1463_fluxblock_Rp02_T6I3.sbml
                for j in range(0,4):
                    file = pref + modelname +'_'+ evaltype + '_R' + str(thr) +'_T'+threshodl+'I'+str(j)+ '.sbml'  # only for inter1
                    if not os.path.exists(file):
                        continue
                    # check if already in overlap
                    if len(overlap_df)>0:
                        currentdf = overlap_df[
                                                (overlap_df['modelname'] == modelname) &
                                                (overlap_df['thr'] == thr) &
                                                (overlap_df['evaltype'] == evaltype) &
                                                (overlap_df['iter'] == j)
                                            ]
                        if currentdf.shape[0]>0:
                            print('already in overlap',flush=True)
                            # if modelname == 'iAF1260' or modelname == 'iNJ661':
                            #     if thr == '2' or thr == '3' or thr == '5' or thr=='p005' or thr=='p001':
                            #         pass
                            #     if j ==1:
                            #         pass
                            # else:
                                # continue
                            continue
                    
                    model = cobra.io.read_sbml_model(file)
                    model_rxns = set(rxn.id.strip('_c') for rxn in model.reactions)
                    # model_mets = set(met.id for met in model.metabolites)
                    # model_genes = set(gene.id for gene in model.genes)
    
                    intersection = len(model_rxns & seedrxns)
                    overlap['modelname'].append(modelname)
                    overlap['thr'].append(thr)
                    overlap['threshod'].append(threshodl)
                    overlap['evaltype'].append(evaltype)
                    overlap['iter'].append(j)
                    overlap['model_subtype'].append('rxns')
                    overlap['both'].append(intersection)
                    overlap['predmodel'].append(len(model_rxns))
                    overlap['biggmodel'].append(len(seedrxns))
                    overlap['ratio_predmodel'].append(intersection/len(model_rxns))
                    overlap['ratio_biggmodel'].append(intersection/len(seedrxns))
                    overlap['tp'].append(intersection)
                    overlap['fp'].append(len(model_rxns)-intersection)
                    overlap['fn'].append(len(seedrxns)-intersection)
                    overlap['precision'].append(intersection/len(model_rxns))
                    overlap['recall'].append(intersection/len(seedrxns))
                    overlap['f1'].append(2*intersection/(len(model_rxns)+len(seedrxns)))
                    

    ## save to csv
    newoverlap_df = pd.DataFrame(overlap)
    print('newoverlap_df:\n',newoverlap_df.shape)
    overlap_df = pd.concat([overlap_df,newoverlap_df],ignore_index=True)
    print(overlap_df.shape)
    # overlap_df.to_csv('model_bigg_overlap_i1.csv',index=False)
    # order by modelname then thr
    overlap_df = overlap_df.sort_values(by=['modelname'])
    # overlap_df.to_csv('model_bigg_overlap_i1.csv',index=False)
    # overlap_df.to_csv('/ibex/user/niuk0a/funcarve/cobra/emodel_bigg_overlap_rerun.csv',index=False)
    # overlap_df.to_csv('/ibex/user/niuk0a/funcarve/cobra/emodel_bigg_overlap_v3.csv',index=False)
    # overlap_df.to_csv(f'/ibex/user/niuk0a/funcarve/cobra/data/result/emodel_bigg_overlap_v4_{evaltype}.csv',index=False)
    overlap_df.to_csv(f'/ibex/user/niuk0a/funcarve/cobra/data/result/emodel_bigg_overlap_v4_{evaltype}_T.csv',index=False)
    return overlap_df
# check_modelbigg(evaltype='allec')
# check_modelbigg(evaltype='flux')

#############################################PLOT###############################################

file_path ='/ibex/user/niuk0a/funcarve/cobra/emodel_bigg_overlap_v4_allec.csv'
df = pd.read_csv(file_path)
## reformat df['thr']
df['thr'] = df['thr'].apply(lambda x: '0.'+x[-1] if x[0]=='p' else x)

# recon_file_path = '/ibex/user/niuk0a/funcarve/cobra/model_bigg_overlap_i1_recon.csv'
recon_data = df[df['evaltype']=='recon']


# data = pd.concat([df, data], ignore_index=True)
data = df
# print(data.head())
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
# evalcol='ratio_biggmodel'   
# evalcolname = 'Ratio of BiGG model reactions' 
# evalcol='f1'   
# evalcolname = 'F1 score' 
evalcol='tp'   
evalcolname = 'True Positive Reactions' 

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
    if evalcol == 'tp':
        refval = int(refval)
        ax.text(1, 0.05, f'Recon={refval}', horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes, color='black', fontsize=10)
    # ax.axhline(y=refval, color='grey', linestyle='--', label='Recon {refval}')
    #set legend for recon reference
    # ax.legend([f'Recon={refval}'], loc='lower left',color='grey')
    else:
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
g.add_legend(title="Rewards")
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
g.set_axis_labels("Iteration", evalcolname)

# 调整整体布局
plt.subplots_adjust(top=0.9)
# g.fig.suptitle("Both vs Iter for Each Modelname", fontsize=16)  # 总标题
# g.fig.suptitle("ratio_biggmodel vs Iter for Each Modelname", fontsize=16)  # 总标题
# g.fig.suptitle("ratio_predmodel vs Iter for Each Modelname", fontsize=16)  # 总标题
g.fig.suptitle(f"Reactions {evalcolname} based on BiGG models", fontsize=16)  # 总标题

# 显示图表
# plt.show()
# plt.savefig('eval_modelrxn_ratio_biggmodel.png', dpi=300)
# plt.savefig('eval_modelrxn_ratio_predmodel.png', dpi=300)
# plt.savefig(f'eval_modelrxn_{evalcol}_recon5.png', dpi=300)
plt.savefig(f'../data/result/figures/eval_modelrxn_{evalcol}_v4_allec.png', dpi=300)
