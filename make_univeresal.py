import cobra
from cobra.io import load_model, load

# 
model = load_model("iYS1720")
# add all bigg reaction to universal model
#curl 'http://bigg.ucsd.edu/api/v2/models'
# curl 'http://bigg.ucsd.edu/api/v2/models/iMM904/genes'
for model in modellist:

