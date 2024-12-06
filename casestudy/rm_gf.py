import cobra
import pandas as pd
gf_m='/Users/niuk0a/Documents/code/gsmm/bvbrc/casestudy/837.83.PATRIC_maxsep_clean_837_gf_tight.sbml'
ngf_m='/Users/niuk0a/Documents/code/gsmm/bvbrc/casestudy/837.83.PATRIC_maxsep_clean_837_ngf_tight.sbml'
ec2r ='/Users/niuk0a/Documents/code/gsmm/recon/reconstructor/reconstructor/Unique_ModelSEED_Reaction_ECs.txt'


ecf = pd.read_csv(ec2r, sep='\t')

gf_model = cobra.io.read_sbml_model(gf_m)
ngf_model = cobra.io.read_sbml_model(ngf_m)
gf_r = gf_model.reactions
ngf_r = ngf_model.reactions

gf_r = set(r.id for r in gf_r)
ngf_r = set(r.id for r in ngf_r)


# rm_rxn = 