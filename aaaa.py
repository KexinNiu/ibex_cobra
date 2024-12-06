f ='/ibex/user/niuk0a/funcarve/cobra/aaseq_CP000148.1_maxsep.csv'
with open(f,'r') as f:
    with open('aaseq_CP000148.1_maxsep_rename.csv','w') as f1:
        count =1
        for line in f:
            line = line.strip('\n')
            a,b = line.split(',',maxsplit=1)
            if len(str(count)) == 1:
                name = 'Gmet_000'+str(count)
            elif len(str(count)) == 2:
                name = 'Gmet_00'+str(count)
            elif len(str(count)) == 3:
                name = 'Gmet_0'+str(count)
            else:
                name = 'Gmet_'+str(count)
            f1.write(name+','+b+'\n')
            count+=1

    

        


