
# from eval_utils import *
# import pandas as pd
# import os
# /ibex/user/niuk0a/CLEAN/app/results/inputs/t2/UP000000268_329726_maxsep.csv


# def get_true_labels(file_name,pref,pred_names,pred_label,pred_probs,pred_all_label):
#     gbkf = pref[:-3] + '.gb'
#     print('gbkf:',gbkf)
#     if not os.path.exists(gbkf):
#         print('no gbkf:!',gbkf)
#         lt2ec = {}
#     # else:
#     #     r_df,df,lt2ec,lt2oldlt, oldlt2lt = read_unif_genebk(file_name,gbkf)

#     r_df = pd.read_csv(file_name,sep='\t')
#     # only take where Reviewed == reviewed
#     r_df = r_df[r_df['Reviewed']=='reviewed']

#     all_label = set()
#     true_label_dict = {}
#     rm_index = []
#     goinfo={}
#     true_label = []
    
#     currpred_names = []
#     currpred_label = []
#     currpred_probs = []
#     print(f'get_true_labels:len(pred_label):{len(pred_label)},len(pred_names):{len(pred_names)},len(pred_all_label):{len(pred_all_label)},len(pred_probs):{len(pred_probs)}')
#     print('pred names:',pred_names[:5])
#     if not lt2ec == {}:
#     #########################################################################
#     #################### WITH GENEBANK + UNIPROT ANNOTATION #################
#         print('WITH GENEBANK + UNIPROT ANNOTATION')
#         for i in range(0,len(pred_names)):
#             # oriname = pred_names[i]
#             # name = oriname.split('|')[1]
#             name = pred_names[i]
#             try:
#                 oldname = lt2oldlt[name]
#                 # print('oldname:',oldname)

#             except:
#                 pass
                
#             remove = True
#             for on in oldname:
#                 ## for each gene names, split by ' ' if there is on , then return unip
#                 # unip  = r_df[r_df['Gene Names'].str.contains(on+' ', na=False) | r_df['Gene Names'].str.endswith(on, na=False)]
#                 unip  = r_df[r_df['Entry'].str.contains(on+' ', na=False)]
#                 # unip = r_df[r_df['Gene Names'].str.contains(on, na=False)]
#                 if on in lt2ec.keys():
#                     if len(lt2ec[name]) > 0:
#                         remove = False
#                         true_ec_lst = lt2ec[name]
#                         for ec in true_ec_lst:
#                             all_label.add(ec)
#                         true_label_dict[name] = true_ec_lst
#                 elif len(unip) != 0:
#                     if unip['EC number'].isnull().values.any():
#                         # print('unip ec number null')
#                         continue
#                     else:
#                         remove = False
#                         uniecs = unip['EC number'].values.tolist()[0]
#                         true_ec_lst = uniecs.split('; ')
#                         for ec in true_ec_lst:
#                             if '-' in ec:
#                                 continue
#                             all_label.add(ec)
#                         true_ec_lst = [ec for ec in true_ec_lst if '-' not in ec]
#                         if len(true_ec_lst) > 0:
#                             true_label_dict[name] = true_ec_lst
#                             # print('name with label found:',name,true_ec_lst)
#                         else:
#                             remove = True
                            
              
#             if remove:
#                 rm_index.append(i)
#             else:
#                 currpred_names.append(name)
#                 currpred_label.append(pred_label[i])
#                 currpred_probs.append(pred_probs[i])
#                 true_label.append(true_label_dict[name])
#     #########################################################################

#     ###################################################################
#     #################### WITH ONLY UNIPROT ANNOTATION #################
#     else:
#         print(' WITH ONLY UNIPROT ANNOTATION ')        
#         for i in range(0,len(pred_names)):
#             # oriname = pred_names[i]
#             # name = oriname.split('|')[1]
#             name = pred_names[i]
#             # unip = r_df[r_df['Gene Names'].str.contains(name, na=False)]
#             unip = r_df[r_df['Entry'].str.contains(name, na=False)]
#             if len(unip) == 0:
#                 # print('no uniprot:',name,uniname)
#                     rm_index.append(i)
#                     continue
#             elif unip['EC number'].isnull().values.any():
#                 rm_index.append(i)
#                 # goinfo[name] = unip['Gene Ontology (GO)']
#                 # continue
#             else:
#                 uniecs = unip['EC number'].values.tolist()[0]
#                 true_ec_lst = uniecs.split('; ')
#                 for ec in true_ec_lst:
#                     all_label.add(ec)
#                 true_label_dict[name] = true_ec_lst
#                 currpred_names.append(name)
#                 currpred_label.append(pred_label[i])
#                 currpred_probs.append(pred_probs[i])
#                 true_label.append(true_label_dict[name])
#     #################### WITH ONLY UNIPROT ANNOTATION #################
#     ###################################################################
#     return true_label, all_label, currpred_names,currpred_label,currpred_probs
       
# def read_unif(f):
#     df = pd.read_csv(f,sep='\t')
#     # only take where Reviewed == reviewed
#     r_df = df[df['Reviewed']=='reviewed'] 
#     #  take columns ['Entry','Protein names','Gene names','Organism','EC number','Status']
#     r_df = r_df[['Entry','EC number']]
#     return r_df,df

# def eval_ori(pref,unif):
#     pred_label, pred_names, pred_all_label = get_pred_labels(pref)
#     pred_names = [name.split('|')[1] for name in pred_names]
#     pred_probs = get_pred_probs(pref)
#     # print('pred names:',len(pred_names),pred_names[:5])
#     true_label, all_label, pred_names, pred_label, pred_probs = get_true_labels(unif,pref,pred_names,pred_label,pred_probs,pred_all_label)
#     return true_label, all_label, pred_names, pred_label, pred_probs




############################GENERATE EVALUATION RESULT#############################
# data_root='/ibex/user/niuk0a/funcarve/cobra/data/tmp'
# resultmatrix = {'modelname':[],'type':[],'iter_round':[],'pre':[],'rec':[],'f1':[],'roc':[],'acc':[],'pred_len':[],'all_len':[]}
# # folder='/ibex/user/niuk0a/funcarve/cobra/data/tmp'   #UP000008561_96561_ORIpredscore_0.pkl'
# for file in os.listdir(data_root):
#     if 'UP' in file and 'ORIpredscore' in file:
#         upid = file.split('_')[0]
#         taxoid = file.split('_')[1]
#         true_label, all_label, pred_names, pred_label, pred_probs = eval_ori(f'/ibex/user/niuk0a/CLEAN/app/results/inputs/t2/{upid}_{taxoid}',f'/ibex/user/niuk0a/funcarve/cobra/uniprot/up_t2/{upid}.tsv')
#         pre,rec,f1,roc,acc = get_eval_metrics(pred_label, pred_probs, true_label, all_label)
#         print('ORI:pre:',pre,'rec:',rec,'f1:',f1,'roc:',roc,'acc:',acc)
#         resultmatrix['modelname'].append(f'{upid}_{taxoid}')
#         resultmatrix['type'].append('CLEAN')
#         resultmatrix['iter_round'].append(0)
#         resultmatrix['pre'].append(pre)
#         resultmatrix['rec'].append(rec)
#         resultmatrix['f1'].append(f1)
#         resultmatrix['roc'].append(roc)
#         resultmatrix['acc'].append(acc)
#         resultmatrix['pred_len'].append(len(pred_names))
#         resultmatrix['all_len'].append(len(all_label))
#         # print(f'{upid}_{taxoid} evaluation done!',resultmatrix)
        
#         orif = f'{data_root}/{upid}_{taxoid}_ORIpredscore_0.pkl'
#         orid = pd.read_pickle(orif)
#         for i in range(1,4):
#             differf = f'{data_root}/{upid}_{taxoid}_differ_predscore_{i}.pkl'
#             if not os.path.exists(differf):
#                 continue
#             # load and evaluate enzyme annotation results
#             differd = pd.read_pickle(differf)
#             currentd = update_differ(orid,differd)

                
#             format_dd = {k.split('|')[1]:v for k,v in currentd.items()}

#             interpred_names,interpred_label_cutoff,interpred_label,interpred_probs = eval_df(format_dd,all_label,pred_names)



#             pre,rec,f1,roc,acc = get_eval_metrics(interpred_label, interpred_probs, true_label, all_label)
#             print(f'IT{i} :pre:',pre,'rec:',rec,'f1:',f1,'roc:',roc,'acc:',acc)
#             resultmatrix['modelname'].append(f'{upid}_{taxoid}')
#             resultmatrix['type'].append('differ')
#             resultmatrix['iter_round'].append(i)
#             resultmatrix['pre'].append(pre)
#             resultmatrix['rec'].append(rec)
#             resultmatrix['f1'].append(f1)
#             resultmatrix['roc'].append(roc)
#             resultmatrix['acc'].append(acc)
#             resultmatrix['pred_len'].append(len(interpred_names))
#             resultmatrix['all_len'].append(len(all_label))
#         # print(f'{upid}_{taxoid} evaluation done!',resultmatrix)

     



# # save evaluation results to csv file   
# file_path = '/ibex/user/niuk0a/funcarve/cobra/data/result/eval_clean_t2UP_specific_ec.csv'

# df = pd.DataFrame(resultmatrix)
# df.to_csv(file_path,index=False)
# print(f'Evaluation results saved to {file_path}')
# exit()


############################PLOT EVALUATION RESULTS#############################
import pandas as pd
import matplotlib.patches as mpatches

f = '/ibex/user/niuk0a/funcarve/cobra/data/result/eval_clean_t2UP_specific_ec.csv'
df = pd.read_csv(f)
clean = df[df['type']=='CLEAN']
differ = df[df['type']=='differ']
## average evaluation results on differ if the have same modelname
differ = differ.groupby('modelname').mean().reset_index()

differ['type'] = 'differ'
## merge clean and differ df
clean = clean.append(differ)
df = clean
## draw plot for different evaluation metrics
import matplotlib.pyplot as plt
import seaborn as sns
## for different 'type' using different color
## draw for precision, recall, f1, roc using box plot
import seaborn as sns
# Initialize a grid of plots with 2 rows and 2 columns
metrics = ["pre", "rec", "f1", "roc"]
fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
colors = sns.color_palette("Set2", len(df["type"].unique()))

# for ax, metric in zip(axes.flat, metrics):
#     sns.boxplot(data=df, x="type", y=metric, ax=ax, palette="Set2")
#     ax.set_title(f"Box Plot of {metric.upper()} by Type")
#     ax.set_xlabel("Type")
#     ax.set_ylabel(metric.upper())
#     xticks = ['CLEAN', 'Ours']
#     ax.set_ylim(0.8, 1)

# draw line plot for each modelname but different evaluation metrics
fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
# colors = sns.color_palette("Set2", len(df["type"].unique()))
colors = sns.color_palette("Set2", len(df["modelname"].unique()))
print('number:',len(df["modelname"].unique()))
for ax, metric in zip(axes.flat, metrics):
    orderdf = df.sort_values(by=["modelname",metric],ascending=False)
    for i, (name, group) in enumerate(orderdf.groupby("type")):
        group = group.sort_values(by=metric,ascending=False)
        # print('group:',group)
        # print(len(group["type"]))
        # ax.plot(group["modelname"], group[metric], label=name, color=colors[i])
        ax.scatter(group["modelname"], group[metric], color=colors[i],marker='d')

    ## add a line to connect the dots for each modelname
    for i, (name, group) in enumerate(orderdf.groupby("modelname")):
        group = group.sort_values(by=metric,ascending=False)
        # if "differ" > "CLEAN" then connect the dots in color "red" else "green"
        if group.iloc[0]["type"] == "differ":
            ax.plot(group["modelname"], group[metric], color=colors[1])
        else:
            ax.plot(group["modelname"], group[metric], color=colors[0])
    
    # ax.legend(labels=["CLEAN", "Ours", "Ours' better than CLEAN", "CLEAN' s better than Ours'"], loc='upper right')

    clean = mpatches.Patch(color=colors[0], label='CLEAN')
    ours = mpatches.Patch(color= colors[1], label='Ours')
    oursb = mpatches.Patch(color=colors[1], label="Ours' better than CLEAN")
    cleanb = mpatches.Patch(color=colors[0], label="CLEAN' s better than Ours'")
    ax.legend(handles=[clean, ours, cleanb,oursb], loc='upper right')


    if metric == "roc":
        mname = "ROC"
    elif metric == "pre":
        mname = "Precision"
    elif metric == "rec":
        mname = "Recall"
    else:
        mname = "F1 Score"
    ax.set_title(f"{mname} by Genomes")
    ax.set_xlabel("Genomes")
    ax.set_ylabel(mname)
    # remove xticks
    ax.set_xticks([])
    ax.set_ylim(min(df[metric])-0.01, max(df[metric])+0.01)





plt.show()
plt.savefig('/ibex/user/niuk0a/funcarve/cobra/data/result/figures/eval_clean_t2UP_specific_ec.png')