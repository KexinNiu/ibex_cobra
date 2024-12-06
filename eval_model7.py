import json
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import cobra
# from matplotlib_venn import venn2
import pickle
prefs =[
    '/ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1',
    '/ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2',
    '/ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2',
    '/ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3_t4',
    '/ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_t4',
    '/ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_t5',
    '/ibex/user/niuk0a/CLEAN/app/results/inputs/AL009126.3_t4'
]
prefs = [i.replace('/ibex/user/niuk0a/CLEAN/app/results/inputs/','/ibex/user/niuk0a/CLEAN/app/results/inputs/sbmls/') for i in prefs]
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
thrlist=['p001','p005','p01','p03','p05','p07','1','2','3','5']

# add on rewards=(0.01 1.5 2.0 3.0 5.0 )
evallist=['allec','flux','fluxblock','block']
# import os
# of = open('eval_model7.err','w')
# thrlist=['p01','p03','p05','p07','1']
# evallist=['allec','flux','fluxblock','block']
# for i in range(0,len(prefs)):
#     pref = prefs[i]
#     modelname = modelnames[i]
#     for thr in thrlist:
#         for evaltype in evallist:
#             for j in range(0,10):
#             # f ='/ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2_maxsep_dfiJN1463_fluxblock_p05iter_10.sbml'
#                 file = pref + '_maxsep_df' + modelname +'_'+ evaltype + '_' + str(thr) +'iter_'+str(j)+ '.sbml'
#                 if not os.path.exists(file):
#                     print(f"File {file} does not exist",file=of)
#                     continue
#                 cmd =f'sbatch memote_sb.sh {file} /ibex/user/niuk0a/CLEAN/app/results/inputs/memote/{modelname}_{evaltype}_{thr}_{i}.json'
#                 os.system(cmd)
            
def check_memotescore(foler):
    #/ibex/user/niuk0a/CLEAN/app/results/inputs/memote

    memotescoredict = {
        'modelname':[],
        'evaltype':[],
        'thr':[],
        'stoichiometric_consistency':[],
        'mass_balance':[],
        'charge_balance':[],
        'metabolite_connectivity':[],
        'unbounded_flux':[],
        'consistency':[]

    }
    print(len(os.listdir(foler)))
    for file in os.listdir(foler):
        if file.endswith('.json'):
            modelname = file.split('_')[0]
            evaltype = file.split('_')[1]
            thr = file.split('_')[2]
            memotescoredict['modelname'].append(modelname)
            memotescoredict['evaltype'].append(evaltype)
            memotescoredict['thr'].append(thr)

            with open(foler+'/'+file) as f:
                data = json.load(f)
                stoichiometric_consistency = 1 - (data['tests']['test_stoichiometric_consistency']['metric'])
                mass_balance = 1 - (data['tests']['test_reaction_mass_balance']['metric'])
                charge_balance = 1 - (data['tests']['test_reaction_charge_balance']['metric'])
                metabolite_connectivity = 1 - (data['tests']['test_find_disconnected']['metric'])
                unbounded_flux = data['tests']['test_find_reactions_unbounded_flux_default_condition']['metric']
    # /ibex/user/niuk0a/CLEAN/app/results/inputs/result_nc.json
                consistency = (stoichiometric_consistency + mass_balance + charge_balance + metabolite_connectivity + (1-unbounded_flux))/5
                memotescoredict['stoichiometric_consistency'].append(stoichiometric_consistency)
                memotescoredict['mass_balance'].append(mass_balance)
                memotescoredict['charge_balance'].append(charge_balance)
                memotescoredict['metabolite_connectivity'].append(metabolite_connectivity)
                memotescoredict['unbounded_flux'].append(unbounded_flux)
                memotescoredict['consistency'].append(consistency)

    df = pd.DataFrame(memotescoredict) 
    ## add refer bigg model information
    # /ibex/user/niuk0a/funcarve/cobra/uniprot/iYO844.xml
    folder='/ibex/user/niuk0a/funcarve/cobra/uniprot'
    for file in os.listdir(folder):
        if file.endswith('.json'):
            modelname = file.split('_')[1]
            evaltype = 'bigg'
            thr = 'bigg'
            memotescoredict['modelname'].append(modelname)
            memotescoredict['evaltype'].append(evaltype)
            memotescoredict['thr'].append(thr)
            with open(folder+'/'+file) as f:
                data = json.load(f)
                stoichiometric_consistency = 1 - (data['tests']['test_stoichiometric_consistency']['metric'])
                mass_balance = 1 - (data['tests']['test_reaction_mass_balance']['metric'])
                charge_balance = 1 - (data['tests']['test_reaction_charge_balance']['metric'])
                metabolite_connectivity = 1 - (data['tests']['test_find_disconnected']['metric'])
                unbounded_flux = data['tests']['test_find_reactions_unbounded_flux_default_condition']['metric']
    biggdf = pd.DataFrame(memotescoredict)
    all_df = pd.concat([df,biggdf])
    all_df = all_df.sort_values(by='modelname')


    #save to csv
    df.to_csv('eval_7results_models.csv',index=False)
    all_df.to_csv('eval_7biggresults.csv',index=False)
    return df
    # return consistency
# df = check_memotescore('/ibex/user/niuk0a/CLEAN/app/results/inputs/memote')
# f='/ibex/user/niuk0a/funcarve/cobra/uniprot/reac_xref.tsv'

def plot_eval(df, svname):
    # 获取所有模型名称
    models = df['modelname'].unique()

    # 设置大图网格
    n_rows = (len(models) + 3) // 4  # 每行最多放4个小图
    fig, axes = plt.subplots(n_rows, 4, figsize=(20, 5 * n_rows), squeeze=False)

    # 遍历每个模型，生成对应的小图
    for idx, model in enumerate(models):
        ax = axes[idx // 4][idx % 4]
        model_data = df[df['modelname'] == model]

        # 获取每种 `evaltype` 数据并绘制柱状图
        categories = ['stoichiometric_consistency', 'mass_balance', 'charge_balance', 
                      'metabolite_connectivity', 'unbounded_flux', 'consistency']
        bar_width = 0.15  # 每条柱的宽度
        x = range(len(categories))  # X 轴为指标位置

        # 遍历每种类型并偏移位置
        for i, (type_name, type_data) in enumerate(model_data.groupby('evaltype')):
            ax.bar(
                [pos + i * bar_width for pos in x], 
                type_data[categories].mean().values,  # 取每种类型的均值
                bar_width, 
                label=type_name
            )

        # 设置子图标题和坐标轴
        ax.set_title(f'Model: {model}')
        ax.set_xticks([pos + bar_width * (len(model_data["evaltype"].unique()) / 2 - 0.5) for pos in x])
        ax.set_xticklabels(categories, rotation=45, ha='right')
        ax.set_ylabel('Score')
        ax.legend()

    # 移除多余的空白子图
    for j in range(len(models), n_rows * 4):
        fig.delaxes(axes[j // 4][j % 4])

    # 调整布局并保存
    plt.tight_layout()
    plt.savefig(svname)
    plt.show()

# 使用示例
# plot_eval(df, 'eval_plotmodel.png')

        
def read_biggmn(file):
    bigg2mn = {}
    mn2seed = {}
        
    with open(file,'r') as f:
        for line in f:
            if line.startswith('*'):
                continue
            else:
                if line.startswith('biggR:'):
                    bigg = line.split('\t')[0].replace('biggR:','')
                    mn = line.split('\t')[1]
                    bigg2mn[bigg] = mn  
                if line.startswith('seedR:'):
                    seed = line.split('\t')[0].replace('seedR:','')
                    mn = line.split('\t')[1]
                    try:
                        mn2seed[mn].append(seed)
                    except:
                        mn2seed[mn] = [seed]
    # save to pickle
        
    with open('bigg2mn.pkl', 'wb') as f:
        pickle.dump(bigg2mn, f)
    with open('mn2seed.pkl', 'wb') as f:
        pickle.dump(mn2seed, f)

    return bigg2mn,mn2seed 

# bigg2mn,mn2seed = read_biggmn('/ibex/user/niuk0a/funcarve/cobra/uniprot/reac_xref.tsv')

             


def check_modelbigg(evaltype='allec'):

    with open('/ibex/user/niuk0a/funcarve/cobra/uniprot/bigg2mn.pkl', 'rb') as f:
        bigg2mn = pickle.load(f)
    with open('/ibex/user/niuk0a/funcarve/cobra/uniprot/mn2seed.pkl', 'rb') as f:
        mn2seed = pickle.load(f)
    outf ='/ibex/user/niuk0a/funcarve/cobra/emodel_bigg_overlap_rerun.csv'
    outf ='/ibex/user/niuk0a/funcarve/cobra/emodel_bigg_overlap_v3.csv'
    if os.path.exists(outf):
        overlap_df = pd.read_csv(outf)
    else:
        overlap_df = pd.DataFrame()
        
    overlap={
        'modelname':[],
        'thr':[],
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
            # 'NC_000853.1_t1_maxsep_dfiLJ478_allec_p01iter_1.sbml'
            print('>>',modelname,thr,flush=True)
            # for j in range(0,10):
            for j in range(0,4):
                file = pref + '_maxsep_df' + modelname +'_'+ evaltype + '_' + str(thr) +'iter_'+str(j)+ '.sbml'  # only for inter1
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
    overlap_df.to_csv('/ibex/user/niuk0a/funcarve/cobra/emodel_bigg_overlap_v3.csv',index=False)
    return overlap_df
check_modelbigg(evaltype='allec')


