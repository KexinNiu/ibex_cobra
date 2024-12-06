name = 'iaf987INTER'
name ='ecoliINTER'
inter = 3
inter= 8
clean_file  ='/ibex/user/niuk0a/funcarve/cobra/aaseq_CP000148.1_maxsep.csv'
clean_file='/ibex/user/niuk0a/funcarve/cobra/NC_000913.3_ecoli_maxsep.csv'

import pickle   


def read_clean_withscore(input_file,threshold=0.8):
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
                    if dis >= threshold:
                        try:
                            pr2ec[pr].append(ecid)
                            # predscore[pr].update({ecid:dis})
                        except KeyError:
                            pr2ec[pr] = [ecid]   
                            # predscore[pr] = {ecid:dis} 
    print('pr2ec-protein number->',len(list(pr2ec.keys())))   
    return pr2ec,ec2pr,predscore

pr2ec,ec2pr,predscore = read_clean_withscore(clean_file,threshold=0.8)
prlist = list(predscore.keys())
allinfo = []
oriinfo = []
for p in prlist:
    info = str(predscore[p])
    oriinfo.append(info)
allinfo.append(oriinfo)
for i in range(1,int(inter)+1):
    flage = [0] * len(prlist)
    nowinfo = oriinfo.copy()
    with open(f'{name}_newpredscore_{i}.pkl', 'rb') as f:
        newpredscore = pickle.load(f)
    with open(f'{name}_updateprs_{i}.pkl', 'rb') as f:
        updateprs = pickle.load(f)
    for p in updateprs:
        flage[prlist.index(p)] = 1
        info = str(newpredscore[p])
        nowinfo[prlist.index(p)] = info
    allinfo.append(nowinfo)
print(len(allinfo))
print(len(allinfo[0]))
print(len(allinfo[1]))
print(len(allinfo[2]))
print(len(prlist))
print(len(allinfo))
# print prlist and allinfo together to a file
with open(f'{name}_allinfo.txt', 'w') as f:
    for i in range(len(prlist)):
        if {allinfo[0][i]}=={allinfo[2][i]}:
            # print(allinfo[0][i],allinfo[-1][i])
            continue
        for j in range (len(allinfo)):
            if j == 0:
                f.write(f'{prlist[i]}\t')
            else:
                f.write(f'{allinfo[j][i]}\t')
        f.write('\n')
        # print('write-->',i)
        # f.write(f'{prlist[i]}\t{allinfo[0][i]}\t{allinfo[1][i]}\t{allinfo[2][i]}\n')
    
        
        

    