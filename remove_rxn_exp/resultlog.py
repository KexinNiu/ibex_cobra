# f ='/ibex/user/niuk0a/funcarve/cobra/result_log.txt'
# f ='/ibex/user/niuk0a/funcarve/cobra/result_log1.txt'
# f = '/ibex/user/niuk0a/funcarve/cobra/result_log_666.txt'
# f = '/ibex/user/niuk0a/funcarve/cobra/result_log_777.txt'
# f = '/ibex/user/niuk0a/funcarve/cobra/result_log_888.txt'
# f = '/ibex/user/niuk0a/funcarve/cobra/result_log_999.txt'
# fs = ['/ibex/user/niuk0a/funcarve/cobra/result_log_555.txt',
#     '/ibex/user/niuk0a/funcarve/cobra/result_log_666.txt',
#     '/ibex/user/niuk0a/funcarve/cobra/result_log_777.txt',
#     '/ibex/user/niuk0a/funcarve/cobra/result_log_888.txt',
#     '/ibex/user/niuk0a/funcarve/cobra/result_log_999.txt'
#     ]
# fs = ['/ibex/user/niuk0a/funcarve/cobra/result_log_555_1.txt',
#     '/ibex/user/niuk0a/funcarve/cobra/result_log_666_1.txt',
#     '/ibex/user/niuk0a/funcarve/cobra/result_log_777_1.txt',
#     '/ibex/user/niuk0a/funcarve/cobra/result_log_888_1.txt',
#     '/ibex/user/niuk0a/funcarve/cobra/result_log_999_1.txt'
#     ]
fs = ['/ibex/user/niuk0a/funcarve/cobra/result_log_555_allunr.txt',
    '/ibex/user/niuk0a/funcarve/cobra/result_log_666_allunr.txt',
    '/ibex/user/niuk0a/funcarve/cobra/result_log_777_allunr.txt',
    '/ibex/user/niuk0a/funcarve/cobra/result_log_888_allunr.txt',
    '/ibex/user/niuk0a/funcarve/cobra/result_log_999_allunr.txt'
    ]
def print_outresult(fs):
    for f in fs:
        # seed = f.split('_')[-1].split('.')[0]
        seed = f.split('_')[2]
        flag = 1
        prec = 0
        print('#######|SEED|',seed,sep='')
        dd =[0,0,0,0,0,0,0,0,0,0,0]
        with open(f, 'r') as f:
            for line in f:
                if line =='\n':
                    continue
                elif line.startswith('######'):
                    continue
                elif line.startswith('# precetage = '):
                    prec = float(line.split(' ')[3])
                    # currentd = [prec]
                elif line.startswith('ori_model FBA'):
                    if flag == 1:
                        # items =['precetage']
                        line = line.strip('\n')
                        titles = line.split('\t')
                        # print('precetage','|',"|".join(titles),sep='')
                        titles.insert(0,'precetage')
                        dd[0] = titles
                        flag = 0
                elif line.startswith('0.04732237605623724'):
                    line = line.strip('\n')
                    items = line.split('\t')
                    # # items = items[:4]+items[5:]
                    # newline = next(f).strip('\n')
                    # newitems = newline.split('\t')
                    # newitems = newitems[:4]+newitems[5:]
                    # print(prec,'|',"|".join(items),sep='')
                    index = int(round(prec/0.05))    
                    # print('index:',index,prec)
                    items.insert(0,str(prec))   
                    dd[index] = items    
                    
        # print (dd)
        for i in dd:
            items = i
            # print('items:',items)
            # try:
            if items != 0:
                print('|'.join(items))
                #
            # print('|'.join(items))

            # except:
            #     pass
            
                
        print('######################')
                    
            
print_outresult(fs)