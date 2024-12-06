
def read_fasta(fasta_file):
    seq=''
    names =[]
    seqs = []
    with open(fasta_file, 'r') as inFile:
        for line in inFile:
            if line.startswith('>'):
                # name = line.strip('\n').split('>')[1]
                name  =line.strip('\n').split('[locus_tag=')[1].split(']')[0]   
                # try:
                #     name = line.strip('\n').split('[gene=')[1].split(']')[0]
                # except IndexError:
                #     name = line.strip('\n').split('[locus_tag=')[1].split(']')[0]   
                names.append(name)
                if seq == '':
                    continue
                else:
                    seqs.append(seq)
                    seq = ''    
            else:
                seq = seq + line.strip('\n')
        seqs.append(seq)
    return names,seqs



f = '/ibex/user/niuk0a/funcarve/cobra/aaseq_CP000148.1.txt'
outf = '/ibex/user/niuk0a/funcarve/cobra/sequence_CP000148.1.fasta'

## only keeps the locus_tag

names, seqs = read_fasta(f)
with open(outf,'w') as f:
    for i in range(len(names)):
        f.write('>'+names[i]+'\n')
        f.write(seqs[i]+'\n')

