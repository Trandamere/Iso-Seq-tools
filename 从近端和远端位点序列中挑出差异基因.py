#coding:utf-8
def FGid(Genelist,Pro,Dis):
    myfile1=open('%s-%s'%(Genelist,Pro),'w')
    myfile2=open('%s-%s'%(Genelist,Dis),'w')
    Gene=[]
    for line1 in open(Genelist):
        Gene.append( line1.strip())
    Title2=[];Fasta2=[]
    for line2 in open(Pro):
        if '>' in line2:
            Title2.append(line2)
        else:
            Fasta2.append(line2)
    for i in Title2:
        if i.split('>')[1].split('-')[0] in Gene:
            print i,Title2.index(i),Fasta2[Title2.index(i)]
            myfile1.write(i)
            myfile1.write(Fasta2[Title2.index(i)])
    Title3=[];Fasta3=[]
    for line3 in open(Dis):
        if '>' in line3:
            Title3.append(line3)
        else:
            Fasta3.append(line3)
    for j in Title3:
        if j.split('>')[1].split('-')[0] in Gene:
            print j,Title3.index(j),Fasta3[Title3.index(j)]
            myfile2.write(j)
            myfile2.write(Fasta3[Title3.index(j)])
    myfile1.close()
    myfile2.close()
    
FGid('12h-0h-down-gene.txt','main-Proxinal-PAS.fasta','main-Distal-PAS.fasta')
