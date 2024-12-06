import pickle
with open('seedr2ec.pkl', 'rb') as f:
    seedr2ec = pickle.load(f)
with open('seedec2r.pkl', 'rb') as f:
    seedec2r = pickle.load(f)
print("seedr2ec dictionary loaded successfully.")
def read_clean_withscore(input_file,threshold=0.5):
    print('threrhold-->',threshold)
    pr2ec = {}
    ec2pr = {}
    predscore = {}
    with open(input_file, 'r') as inFile:
        for line in inFile:
            line = line.strip('\n')
            line = line.split(',')
            pr = line[0]
            # items = line[-1].split(',')
            items = line[1:]
            for item in items:
                if item.startswith('EC:'):
                    ec,dis = item.split('/')
                    ecid = ec.split(':')[-1]
                    # ecid = ec
                    dis = float(dis)
                    if dis >= 0.0001:
                        try:
                            predscore[pr].update({ecid:dis})
                        except:
                            predscore[pr] = {ecid:dis}
                        try:
                            ec2pr[ecid].update({pr:dis})
                        except:
                            ec2pr[ecid] = {pr:dis}
                    # if dis >= threshold:
                    if dis <= threshold: ### reverse the threshold
                        try:
                            pr2ec[pr].append(ecid)
                            # predscore[pr].update({ecid:dis})
                        except KeyError:
                            pr2ec[pr] = [ecid]   
                            # predscore[pr] = {ecid:dis} 
    print('pr2ec-protein number->',len(list(pr2ec.keys())))   
    return pr2ec,ec2pr,predscore

def update_predscore(newreactions,ec2pr,newpredscore,updateprs,reward):
    for r in newreactions:
        r = r.split('_')[0]
        try:
            ecid = seedr2ec[r]
            print('ecid:', ecid)
        except KeyError:
            print('no ecid:next')
            continue
        for ec in ecid:
            try:
                prd = ec2pr[ec]
                print('prd:', prd)
            except KeyError:
                print('no prd:next')
                continue

            # take the pr with max val
            try:
                pr = max(prd, key=prd.get)
                updateprs.append(pr)
                print('pr:', pr)
            except ValueError:
                print('no pr:next')
                continue

            try:
                # s = min(float(newpredscore[pr][ec]) + reward, 1)
                s = max(float(newpredscore[pr][ec]) + reward, 0)
                print('s:', s)
                newpredscore[pr].update({ec: s})
                print('add+')
            except KeyError:
                print('no score:next')
                    
    return newpredscore,updateprs


clean_file = '/ibex/user/niuk0a/funcarve/cobra/aaseq_CP000148.1_maxsep.csv'        
pr2ec,ec2pr,predscore= read_clean_withscore(clean_file,threshold=0.5)
# print('pr2ec:',pr2ec) #pr2ec: {'gyrB': ['EC:5.6.2.2'], 'gyrA': ['EC:5.6.2.2'], 'gpsA': ['EC:1.1.1.94'], 'fgrL': ['EC:2.7.13.3'], 'Gmet_0015': ['EC:1.3.98.6'], 'hemE': ['EC:4.1.1.37'], 'hemH': ['EC:4.99.1.1'], 'groEL': ['EC:5.6.1.7'], 'Gmet_0042': ['EC:3.5.1.9'], 'cysS': ['EC:6.1.1.16'], 'glnS': ['EC:6.1.1.18'], 'ispF': ['EC:4.6.1.12'], 'ispD': ['EC:2.7.7.60'], 'selA': ['EC:2.9.1.1'], 'rpe': ['EC:5.1.3.1'], 'Gmet_0069': ['EC:2.7.7.65'], 'glnE': ['EC:2.7.7.89', 'EC:2.7.7.42'], 'mtnA': ['EC:5.3.1.23'], 'gatA': ['EC:6.3.5.7'], 'Gmet_0077': ['EC:2.7.13.3'], 'Gmet_0085': ['EC:3.1.3.48'], 'ogt': ['EC:2.1.1.63'], 'glmU': ['EC:2.7.7.23', 'EC:2.3.1.157'], 'glmS-1': ['EC:2.6.1.16'], 'hypE': ['EC:2.7.9.3'], 'pyk': ['EC:2.7.1.40'], 'ppk': ['EC:2.7.4.1'], 'Gmet_0135': ['EC:5.4.2.2'], 'ppiA': ['EC:5.2.1.8'], 'nuoK-2': ['EC:1.12.98.3'], 'nuoM-2': ['EC:1.12.98.3'], 'tkt': ['EC:2.2.1.1'], 'Gmet_0188': ['EC:3.1.26.11'], 'hemB': ['EC:4.2.1.24'], 'Gmet_0194': ['EC:4.99.1.12'], 'Gmet_0197': ['EC:2.7.13.3'], 'alaS': ['EC:6.1.1.7'], 'argB': ['EC:2.7.2.8'], 'argD': ['EC:2.6.1.11'], 'argF': ['EC:2.1.3.3'], 'argG': ['EC:6.3.4.5'], 'argH': ['EC:4.3.2.1'], 'dapA': ['EC:4.3.3.7'], 'dapB': ['EC:1.17.1.8'], 'dapL': ['EC:2.6.1.83'], 'folD-1': ['EC:3.5.4.9', 'EC:1.5.1.5'], 'Gmet_0256': ['EC:6.2.1.44'], 'cysM': ['EC:2.5.1.47'], 'truC': ['EC:5.4.99.26'], 'msrA': ['EC:1.8.4.11'], 'ltaA': ['EC:4.1.2.48', 'EC:4.1.2.49'], 'Gmet_0300': ['EC:2.7.7.77'], 'moaA-1': ['EC:4.1.
# exit()
# print('ec2pr:',ec2pr) #'EC:1.6.5.7': {'Gmet_3556': 0.0001}, 'EC:2.1.1.17': {'Gmet_3557': 0.0007}}
# exit()
# print('predscore:',predscore) #8': 0.0001, 'EC:3.2.2.31': 0.0001}, 'Gmet_3542': {'EC:3.1.11.6': 0.0062}, 'nadA': {'EC:2.5.1.72': 0.9919}, 'yrdA': {'EC:2.3.1.197': 0.1841}, 'Gmet_3546': {'EC:1.16.1.7': 0.0002}, 'mfd': {'EC:3.6.4.12': 0.0003}, 'Gmet_3548': {'EC:5.2.1.8': 0.7986}, 'Gmet_3549': {'EC:5.2.1.8': 0.9802}, 'Gmet_3550': {'EC:5.2.1.8': 0.004}, 'hemY': {'EC:1.3.3.15': 0.9626}, 'Gmet_3552': {'EC:1.14.99.60': 0.0013}, 'def-2': {'EC:3.5.1.88': 0.9892}, 'Gmet_3554': {'EC:2.2.1.6': 0.0146}, 'Gmet_3556': {'EC:1.6.5.7': 0.0001}, 'Gmet_3557': {'EC:2.1.1.17': 0.0007}, 'gidB': {'EC:2.1.1.170': 0.9862}, 'mnmE': {'EC:1.2.1.70': 0.0001}, 'Gmet_3562': {'EC:1.12.98.3': 0.0001}, 'rnpA': {'EC:3.1.26.5': 0.9919}}
# exit()
new_reactions =  {'rxn05241_c', 'rxn01396_c', 'rxn02297_c', 'rxn01465_c', 'rxn34742_c', 'rxn02832_c', 'rxn13893_c', 'rxn07466_c', 'rxn40643_c', 'rxn03880_c', 'lipid_rxn', 'rxn20288_c', 'rxn03167_c', 'rxn26474_c', 'biomass_GmPos', 'rxn41888_c', 'rxn09106_c', 'rxn08240_c', 'rxn12787_c', 'rxn10280_c', 'protein_rxn', 'dna_rxn', 'rxn03866_c', 'rxn41785_c', 'rxn02010_c', 'rxn01790_c', 'rxn05468_c', 'rxn03030_c', 'rxn32502_c', 'rxn35484_c', 'rxn19742_c', 'rxn39267_c', 'rxn01972_c', 'rxn08801_c', 'rxn44825_c', 'rxn30741_c', 'rxn05301_c', 'rxn16332_c', 'rxn23730_c', 'rxn13699_c', 'rxn31476_c', 'rxn19840_c', 'rxn40439_c', 'rxn13697_c', 'rxn08797_c', 'rxn46214_c', 'rxn00979_c', 'EX_biomass', 'rxn42058_c', 'rxn09631_c', 'rxn03617_c', 'rxn02866_c', 'teichoicacid_rxn', 'rxn05039_c', 'rxn07189_c', 'rxn08186_c', 'rxn00940_c', 'rxn13698_c', 'rxn26011_c', 'rxn09639_c', 'rxn09865_c', 'rxn03108_c', 'rxn01476_c', 'rxn26471_c', 'rxn00777_c', 'rxn31310_c', 'rxn00621_c', 'rxn40487_c', 'rxn01258_c', 'rxn26475_c', 'rxn08019_c', 'rxn10281_c', 'rxn01504_c', 'peptidoglycan_rxn', 'rxn00963_c', 'rxn34741_c', 'rxn39069_c', 'rna_rxn', 'rxn12249_c', 'rxn08921_c', 'rxn08086_c', 'rxn42988_c', 'rxn00346_c', 'rxn08856_c', 'rxn05144_c', 'rxn03086_c', 'rxn34424_c', 'rxn05300_c', 'rxn08799_c', 'rxn01538_c', 'rxn10272_c', 'rxn18910_c', 'rxn24330_c', 'rxn10626_c', 'rxn42899_c', 'EX_cpd00559_e', 'rxn18920_c', 'rxn05663_c', 'rxn20769_c', 'rxn47859_c', 'rxn16597_c', 'rxn30323_c', 'rxn47471_c', 'rxn10271_c', 'rxn19774_c', 'rxn08088_c', 'rxn05226_c', 'rxn02898_c', 'rxn37182_c', 'rxn00340_c', 'rxn00122_c', 'rxn02304_c', 'rxn05220_c', 'rxn14012_c', 'rxn40745_c', 'rxn27379_c', 'rxn08233_c', 'GmPos_cellwall', 'rxn03874_c', 'rxn20782_c', 'rxn47074_c', 'rxn38026_c', 'rxn05242_c', 'cofactor_rxn', 'rxn05618_c', 'rxn03084_c', 'rxn30368_c', 'rxn08084_c', 'rxn03852_c'}
# print('seedr2ec:',seedr2ec)
# exit()
newpredscore = predscore.copy()
updateprs = []
reward = 0.1
newpredscore, updateprs = update_predscore(new_reactions,ec2pr,newpredscore,updateprs,reward)
print('updateprs:',updateprs)
