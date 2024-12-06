import pickle
import pandas as pd


def read_pickle(f):
    with open(f, 'rb') as handle:
        data = pickle.load(handle)
    return data

data = read_pickle('/ibex/user/niuk0a/funcarve/cobra/iHN637_R_newpredscore_1.pkl')
data9 = read_pickle('/ibex/user/niuk0a/funcarve/cobra/iHN637_R_newpredscore_5.pkl')

for key, item in data.items():
    if key in data9:
        if item != data9[key]:
            print(key, item, data9[key])
    else:
        continue