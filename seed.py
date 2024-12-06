f = '/ibex/user/niuk0a/funcarve/cobra/SEED2VMH_translation.csv'

def make_dict(f):
    seed2bigg={}
    bigg2seed={}
    with open(f, 'r') as f:
        for line in f:
            line = line.strip().split(',')
            seed = line[0]
            bigg = line[1]
            if seed.startswith('EX_'):
                seed = seed.split('_')[1]
                seed = seed.split('(')[0]
                seed = seed+'_e'
                if bigg.startswith('EX_'):
                    bigg = bigg.split('_')[1]
                    bigg = bigg.split('(')[0]
                    bigg = bigg+'_e'
            seed2bigg[seed] = bigg
            bigg2seed[bigg] = seed
    return seed2bigg,bigg2seed

seed2bigg,bigg2seed = make_dict(f)
## save 
import pickle
with open('seed2bigg.pkl', 'wb') as f:
    pickle.dump(seed2bigg, f)
with open('bigg2seed.pkl', 'wb') as f:
    pickle.dump(bigg2seed, f)
            
            
            