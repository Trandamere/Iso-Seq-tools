#coding:utf-8
import os
def lncRNA(gff):
    '''
    #1 与FG41比较找出新的转录本
    os.system(r'/home/luping/gffcompare-0.10.2/gffcompare -r FGRAMPH1-41version-cdna.gff3 -o ./FG41-%s %s'%(gff,gff))
    #2 去除SJ相同或，在编码基因区域的转录本= c e j k

    SameSJ=['=','c','e','j','k']
    NewT=[]
    for line1 in open('FG41-%s.tracking'%gff):
        if line1.split('\t')[3] not in SameSJ:
            T= line1.split('\t')[4].split('|')[1] 
            NewT.append(T)
    #3 提取NewT转录本
    myfile=open('NewT-%s'%gff,'w')
    for Tid in NewT:
        print Tid
        for line2 in open(gff):
            if '%s"'%Tid in line2:
                myfile.write(line2)
    myfile.close()

    #4 提取序列进行PLEK预测
    
    os.system('gffread -w %s.fasta -g ph1.fasta NewT-%s'%(gff,gff))

    #5 进行PLEK预测
    os.system('python /disk/luping/bam/iso-bam/PLEK.1.2/PLEK.py -fasta %s.fasta -out lnc-%s.fasta -thread 1'%(gff,gff))
    '''
    #6 去除lncRNA序列进行进行Pfam和Rfam预测
    PLEK=[]
    lncFile='lnc-%s.fasta'%gff
    for line3 in open(lncFile):
        if 'Non-coding' in line3:
            PLEK.append( line3.split('>')[1].split(' ')[0])
    #选出gff提取序列
    myfile1=open('PLEK-%s'%gff,'w')
    for PLEKid in PLEK:
        print PLEKid
        for line4 in open(gff):
            if '%s"'%PLEKid in line4:
                myfile1.write(line4)
    myfile1.close()
    #提取序列
    os.system('gffread -w PLEK-%s.fasta -g ph1.fasta PLEK-%s'%(gff,gff))
    
#lncRNA('Fusion-PH1-fusion-version2.3.gtf.gff3.gtf')
lncRNA('Fusion-embossORF-gene-delclassS-C-C-FilterPT-Hy-fusion.gff.gff3.gtf')
lncRNA('Fusion-embossORF-gene-delclassS-C-C-FilterPT-Sex-fusion.gff.gff3.gtf')
