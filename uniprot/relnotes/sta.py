'''
UniProtKB Release 2011_01 consists of 13,069,501 entries (UniProtKB/Swiss-Prot: 
524,420 entries and UniProtKB/TrEMBL: 12,545,081 entries)
UniRef100 Release 2011_01 consists of 11,659,891 entries
UniRef90 Release 2011_01 consists of 7,623,063 entries
UniRef50 Release 2011_01 consists of 3,653,743 entries
UniParc Release 2011_01 consists of 24,828,830 entries
UniMES Release 1.0 consists of 6,028,191 entries
'''
def readnote(fp):
    with open(fp,'r') as f:
        lines = f.readlines()
        entries = int(str(lines[8].split()[5]).replace(',',''))
        swissp = int(str(lines[9].split()[0]).replace(',',''))
        trembl = int(str(lines[9].split()[-2]).replace(',',''))
        ref100 = int(str(lines[10].split()[-2]).replace(',',''))
        ref90 = int(str(lines[11].split()[-2]).replace(',',''))
        ref50 = int(str(lines[12].split()[-2]).replace(',',''))
        uparc = int(str(lines[13].split()[-2]).replace(',',''))
    return  entries,swissp,trembl,ref100,ref90,ref50, uparc
    

for year in range(11, 25):
    # print('year:',year)
    f = f'/ibex/user/niuk0a/funcarve/cobra/uniprot/relnotes/uniprot_{year}_01.txt'
    # print(f)
    entries,swissp,trembl,ref100,ref90,ref50,uparc = readnote(f)
    year = '20'+str(year)
    print(year,entries,swissp,trembl,ref100,ref90,ref50, uparc,sep='|')
    # 2024 uprac changed to 607,912,929
    
    