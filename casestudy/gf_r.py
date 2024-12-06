import cobra
import pandas as pd
gf_m='/ibex/user/niuk0a/funcarve/cobra/casestudy/837.83.PATRIC_maxsep_clean_837_gf_tight.sbml'
ngf_m='/ibex/user/niuk0a/funcarve/cobra/casestudy/casestudy/837.83.PATRIC_maxsep_clean_837_ngf_tight.sbml'
ec2r ='/ibex/user/niuk0a/funcarve/reconstructor/reconstructor/Unique_ModelSEED_Reaction_ECs.txt'

ecf = pd.read_csv(ec2r, sep='\t')

gf_model = cobra.io.read_sbml_model(gf_m)
ngf_model = cobra.io.read_sbml_model(ngf_m)
gf_r = gf_model.reactions
ngf_r = ngf_model.reactions

gf_r = set(r.id for r in gf_r)
ngf_r = set(r.id for r in ngf_r)

new_r = gf_r - ngf_r
print('newr_number:',len(new_r)) # 186 
new_r2ec ={}
notfound = []
for r in new_r:
    r = r.replace('_c','')
    try:
        ec = ecf[ecf['ModelSEED ID'] == r]['External ID'].values
        # print(r,ec)
        new_r2ec[r] = ec[0]
    except:
        notfound.append(r)
        continue
print(new_r2ec) 
print(len(new_r2ec)) # 69 
#### 'rxn01603': '6.3.2.12', 刚好below threshold
#fig|837.83.peg.1594 Dihydrofolate synthase (EC 6.3.2.12) @ Folylpolyglutamate synthase (EC 6.3.2.17),EC:6.3.2.12/0.7038

# {'rxn09111': '2.7.8.5', 'rxn24055': '4.2.3.52', 'rxn05618': '2.A.55.-.-', 'rxn26435': '4.1.2.-', 'rxn05040': '4.1.99.12', 
# 'rxn10204': '2.3.1.15', 'rxn00774': '6.3.4.7', 'rxn15287': '2.7.1.1', 'rxn25213': '2.4.1.-', 'rxn33860': '5.3.1.6', 
# 'rxn02937': '6.3.3.1', 'rxn46500': '1.14.13.M33', 'rxn01396': '2.7.1.-', 'rxn15161': '5.1.3.8', 'rxn00122': '2.7.7.2', 
# 'rxn39451': '1.17.1.8', 'rxn03086': '2.6.1.-', 'rxn01603': '6.3.2.12', 'rxn37219': '2.3.1.51', 'rxn35957': '1.1.1.71', 
# 'rxn23043': '2.5.1.3', 'rxn43849': '2.6.1.111', 'rxn10202': '2.3.1.15', 'rxn00293': '2.7.7.23', 'rxn10212': '2.3.1.51', 
# 'rxn36480': '1.2.1.5', 'rxn19277': '1.17.4.1', 'rxn03866': '2.3.1.-', 'rxn42899': '3.5.4.-', 'rxn03030': '2.3.1.89', 
# 'rxn01802': '1.3.99.7', 'rxn15158': '2.3.1.3', 'rxn09113': '2.7.8.5', 'rxn03141': '6.3.2.10', 'rxn03108': '2.7.4.7', 
# 'rxn00137': '3.1.3.7', 'rxn01974': '5.1.1.7', 'rxn02010': '6.3.2.7', 'rxn02130': '3.5.1.92', 'rxn00309': '1.4.1.18', 
# 'rxn01439': '2.7.1.1', 'rxn09106': '3.1.3.27', 'rxn01467': '4.2.3.10', 'rxn00621': '2.7.7.39', 'rxn10211': '2.3.1.51', 
# 'rxn00891': '3.5.1.33', 'rxn38902': '3.6.1.15', 'rxn00616': '1.1.5.3', 'rxn34055': '6.3.5.3', 'rxn26436': '2.7.1.-', 
# 'rxn18920': '1.14.13.27', 'rxn10203': '2.3.1.15', 'rxn05663': '2.A.3.1.-', 'rxn19774': '2.7.8.33', 'rxn03617': '1.1.1.259', 
# 'rxn44067': '2.4.1.du', 'rxn01465': '3.5.2.3', 'rxn07295': '1.4.3.19', 'rxn01972': '3.5.1.47', 'rxn41785': '1.17.1.11', 
# 'rxn35438': '3.6.1.41', 'rxn01513': '2.7.4.12', 'rxn40439': '2.6.1.105', 'rxn03880': '4.2.1.-', 'rxn00461': '2.5.1.7', 
# 'rxn18910': '3.6.3.9', 'rxn15164': '3.2.2.14', 'rxn05301': '2.A.3.1.-', 'rxn09109': '2.7.8.5'}

print('notfound:',notfound) # 117
# notfound: ['EXpd00051_e', 'rxn32546', 'EXpd00027_e', 'EXpd00066_e', 'EXpd00129_e', 'EXpd00099_e', 'rxn25583', 'EXpd01160_e', 'rxn05242', 'rxn08188', 'rxn05241', 'rxn46587', 'EXpd28082_e', 'EXpd00035_e', 'rxn05221', 'EXpd01080_e', 'rxn12249', 'rxn41088', 'EXpd00039_e', 'rxn10272', 'rxn09632', 'EXpd15237_e', 'EXpd00132_e', 'rxn30054', 'rxn34357', 'EXpd00054_e', 'rxn30508', 'teichoicacid_rxn', 'EXpd00322_e', 'rxn46266', 'rxn10271', 'rxn09306', 'EXpd00001_e', 'rxn10281', 'EXpd15269_e', 'rxn09467', 'EXpd00084_e', 'rxn12622', 'rxn42058', 'rxn28184', 'EXpd00067_e', 'cofactor_rxn', 'EXpd00254_e', 'EXpd03847_e', 'rxn39245', 'rxn05573', 'EXpd00053_e', 'rxn08921', 'EXpd00161_e', 'dna_rxn', 'EXpd00834_e', 'EXpd01107_e', 'EXpd00033_e', 'EXpd00205_e', 'rxn27379', 'EXpd00149_e', 'rxn09934', 'rxn42988', 'rxn05226', 'EXpd00069_e', 'rxn08131', 'EXpd00070_e', 'EXpd00010_e', 'rxn28181', 'rxn47859', 'EXpd03846_e', 'rxn10030', 'biomass_GmNeg', 'EXpd00107_e', 'rxn40745', 'rxn05305', 'rxn08233', 'rxn23050', 'rxn29048', 'EXpd28301_e', 'EXpd00063_e', 'rxn25238', 'rxn09532', 'rxn05216', 'rxn44825', 'rna_rxn', 'EXpd00119_e', 'EXpd00060_e', 'EXpd00023_e', 'rxn08240', 'EXpd00009_e', 'EXpd00214_e', 'rxn10878', 'EXpd00041_e', 'rxn12385', 'rxn46214', 'rxn32282', 'rxn47813', 'rxn47074', 'lipid_rxn', 'EXpd00030_e', 'rxn39472', 'EX_biomass', 'EXpd00156_e', 'EXpd01113_e', 'rxn12854', 'rxn34491', 'rxn08982', 'EXpd00971_e', 'EXpd10515_e', 'EXpd00075_e', 'rxn47471', 'peptidoglycan_rxn', 'EXpd17026_e', 'EXpd00065_e', 'rxn28180', 'rxn28183', 'protein_rxn', 'EXpd01741_e', 'rxn10643', 'rxn10280', 'EXpd15298_e']

clean_pref = '/Users/niuk0a/Documents/code/gsmm/bvbrc/casestudy/837.83.PATRIC_maxsep.csv'
predscore={}
def read_clean(input_file):
    threrhold = 0
    print('threrhold-->',threrhold)
    pr2ec = {}
    with open(input_file, 'r') as inFile:
        for line in inFile:
            line = line.split()
            pr = line[0]
            items = line[-1].split(',')
            for item in items:
                if item.startswith('EC:'):
                    ec,dis = item.split('/')
                    ecid = ec.split(':')[-1]
                    dis= float(dis)
                    if dis >= threrhold:
                        try:
                            pr2ec[pr].append(ecid)
                            predscore[pr].update({ecid:dis})
                        except KeyError:
                            pr2ec[pr] = [ecid]   
                            predscore[pr] = {ecid:dis} 
    print('pr2ec-protein number->',len(list(pr2ec.keys())))   
    return pr2ec,predscore
pr2ec ,predscore = read_clean(clean_pref)
# make ec2pr and only keeps the high confidence protein and recode the score
ec2pr = {}
scored = {}
for pr,ecs in pr2ec.items():
    for ec in ecs:
        if ec in ec2pr.keys():
            if predscore[pr][ec] > scored[ec]:
                ec2pr[ec] = [pr]
        else:
            ec2pr[ec] = [pr]
            scored[ec] = predscore[pr][ec]
total = 0
justb = 0
for rxn in new_r2ec.keys():
    ec = new_r2ec[rxn]
    if ec in ec2pr.keys():
        # print('if pr in predscore:',)
        total +=1
        if scored[ec] > 0.6:
            justb +=1
            print('correct:',rxn,ec,'>',ec2pr[ec],'|',scored[ec])
            # print(rxn,ec,'>',ec2pr[ec],'|',scored[ec])
        # print(rxn,ec,'>',ec2pr[ec],'|',scored[ec])
    else:
        total +=1 # total 69 reaction 
        # print(rxn,ec,'notfound')

        ## among the gapfilling reactions, 
        # 4 of them have a protein with score > 0.6, but below threshold 0.8
print('total:',total)#total: 24/69(should have a enzyme/protein mapping to this rxn in the model)
print ('justb:',justb)#justb: 4 (there is a protein with score > 0.6, but below threshold 0.8)
print('ratio:',justb/total)# ratio: 0.16666666666666666 




# correct: rxn03141 6.3.2.10 > ['fig|837.83.peg.952'] | 0.7706
# correct: rxn01603 6.3.2.12 > ['fig|837.83.peg.1594'] | 0.7038 
# Methyl-coenzyme M reductase I, protein D 2.8.4.1 
# ATP + L-Glutamate + Dihydropteroate => ADP + Phosphate + H+ + Dihydrofolate
# correct: rxn00461 2.5.1.7 > ['fig|837.83.peg.632'] | 0.6644
# correct: rxn00137 3.1.3.7 > ['"fig|837.83.peg.966'] | 0.7893
# H2O + Adenosine 3-5-bisphosphate <=> Phosphate + AMP