gff='/ibex/user/niuk0a/funcarve/bvbrc/memote/kegg_stat_gf.csv'
ngff='/ibex/user/niuk0a/funcarve/bvbrc/memote/kegg_stat_ngf.csv'


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def compare(gff,ngff):
    gft = pd.read_csv(gff)
    ngft = pd.read_csv(ngff)
    print('gft:',gft.shape)
    print('ngft:',ngft.shape)
    # keeps only the common taxo_id
    common_taxo_id = pd.merge(gft,ngft,on='taxo_id',how='inner')
    print('common_taxo_id:',common_taxo_id.shape)
    print(common_taxo_id.columns)
    print(common_taxo_id.head())
    totalx = common_taxo_id['keggout_total_score_x']
    totaly = common_taxo_id['keggout_total_score_y']
    total_met_x = common_taxo_id['keggout_total_metabolites_x']
    total_met_y = common_taxo_id['keggout_total_metabolites_y']
    total_rxn_x = common_taxo_id['keggout_total_reactions_x']
    total_rxn_y = common_taxo_id['keggout_total_reactions_y']
    total_genes_x = common_taxo_id['keggout_total_genes_x']
    total_genes_y = common_taxo_id['keggout_total_genes_y']
    total_compartments_x = common_taxo_id['keggout_total_compartments_x']
    total_compartments_y = common_taxo_id['keggout_total_compartments_y']
    metabolic_coverage_x = common_taxo_id['keggout_metabolic_coverage_x']
    metabolic_coverage_y = common_taxo_id['keggout_metabolic_coverage_y']
    unconserved_metabolites_x = common_taxo_id['keggout_unconserved_metabolites_x']
    unconserved_metabolites_y = common_taxo_id['keggout_unconserved_metabolites_y']
    consistency_x = common_taxo_id['keggout_consistency_x']
    consistency_y = common_taxo_id['keggout_consistency_y']
    annotation_met_x = common_taxo_id['keggout_annotation_met_x']
    annotation_met_y = common_taxo_id['keggout_annotation_met_y']
    annotation_rxn_x = common_taxo_id['keggout_annotation_rxn_x']
    annotation_rxn_y = common_taxo_id['keggout_annotation_rxn_y']
    annotation_gene_x = common_taxo_id['keggout_annotation_gene_x']
    annotation_gene_y = common_taxo_id['keggout_annotation_gene_y']
    annotation_sbo_x = common_taxo_id['keggout_annotation_sbo_x']
    annotation_sbo_y = common_taxo_id['keggout_annotation_sbo_y']
    
    plt.figure(figsize=(10, 5))
    plt.boxplot([totalx,totaly],labels=['Gapfill', 'No Gapfill'])
    plt.title('Boxplot of Total Score')
    plt.ylabel('Total Score')
    # plt.show()
    plt.savefig('kegg_total_score.png')
    
    # #new plt
    plt.figure(figsize=(10, 5))
    # plt.subplot(1, 2, 1)
    plt.boxplot([totalx,totaly,
                #  total_compartments_x,total_compartments_y,total_genes_x,total_genes_y,total_met_x,total_met_y,total_rxn_x,total_rxn_y,metabolic_coverage_x,metabolic_coverage_y,
                #  unconserved_metabolites_x,unconserved_metabolites_y,
                 consistency_x,consistency_y,annotation_met_x,annotation_met_y,annotation_rxn_x,annotation_rxn_y,annotation_gene_x,annotation_gene_y,annotation_sbo_x,annotation_sbo_y],
                labels=['Total Score_gf','_ngf', 
                        # 'Total Compartments_gf','Total Compartments_ngf','Total Genes_gf','Total Genes_ngf','Total Metabolites_gf','Total Metabolites_ngf','Total Reactions_gf',
                        # 'Total Reactions_ngf','Metabolic Coverage_gf','Metabolic Coverage_ngf','Unconserved Metabolites_gf','Unconserved Metabolites_ngf',
                        'Consistency_gf','_ngf','Met_gf','_ngf','Rxn_gf','_ngf','Gene_gf','_ngf','Sbo_gf','_ngf'])
                        
    plt.title('Boxplot of KEGGout Results')
    plt.ylabel('Scores')
    # plt.show() 
    plt.savefig('kegg_gf_ngf.png')
    
    #new plt
    plt.figure(figsize=(10, 5))
    plt.boxplot(
        [total_compartments_x,total_compartments_y,total_genes_x,total_genes_y,total_met_x,total_met_y,total_rxn_x,total_rxn_y,metabolic_coverage_x,metabolic_coverage_y,
                 unconserved_metabolites_x,unconserved_metabolites_y],
        labels=[
            # 'Total Compartments_gf','Total Compartments_ngf','Total Genes_gf','Total Genes_ngf','Total Metabolites_gf','Total Metabolites_ngf','Total Reactions_gf',
            #             'Total Reactions_ngf','Metabolic Coverage_gf','Metabolic Coverage_ngf','Unconserved Metabolites_gf','Unconserved Metabolites_ngf',
            'Compart_gf','_ngf','Genes_gf','_ngf','Met_gf','_ngf','Rxn_gf','_ngf','MetaCoverage_gf','_ngf',
            'UnconservedMeta_gf','_ngf']
    )
    plt.title('Boxplot of KEGGout_met_rxn Results')
    plt.ylabel('Number of Metabolites/Reactions')
    # plt.show()
    plt.savefig('kegg_met_rxn.png')
    return
# plt.boxplot([totalx, totaly], labels=['Gapfill', 'No Gapfill'])
compare(gff,ngff)