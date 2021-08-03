#coding:utf-8
#luping
#2019.2.22
#进行转录本T0排序
import os
#将PB,RNAseq与FG的注释分离。PB，RNAseq进行ORF预测，FG的用原注释

def twoPart(x):
    myfile1=open(r'PB-%s'%x,'w')
    myfile2=open(r'FG-%s'%x,'w')
    for line1 in open(r'%s'%x):
        if 'PH1\t' in line1:
            myfile2.write(line1)
        else:
            myfile1.write(line1)
    myfile1.close()
    myfile2.close()
#twoPart('replaceID-ReFG-StringTieNew-StingTie1To1-SuppShort-add-Del-classcodeS-delEXP-delLess200bp-C-C-FilterPT-ISO-fusion.collapsed.gff')
#twoPart('replaceID-ReFG-StringTieNew-StingTie1To1-SuppShort-add-Del-classcodeS-delEXP-delLess200bp-C-C-FilterPT-ISO-nofusion.collapsed.gff')
#使用transcript-T0进行PB ORF预测
import os
def joinfasta(x,y):#连接fasta序列
    myfile=open(r'%s'%y,'w')
    for line in open(r'%s'%x):
        if '>' in line:
            #print line
            myfile.write('\n%s'%line)
        else:
            #print line.strip()
            myfile.write(line.strip())
    myfile.close()
def join2fasta(x,y):#去除第一行空白
    myfile=open(r'%s'%y,'w')
    for line in open(r'%s'%x):
        if line=='\n':
            pass
        else:
            myfile.write(line)
    myfile.close()
def derev(x,y):#去除反向ORF
    myfile=open(r'%s'%y,'w')
    title=[]
    fasta=[]
    for line in open(r'%s'%x):
        if '>' in line:
            title.append(line)
        else:
            fasta.append(line)
    while len(title)>0:
        titlepop=title.pop()
        fastapop=fasta.pop()
        if 'REVERSE SENSE' in titlepop:
            pass
        else:
            myfile.write(titlepop)
            myfile.write('%s\n'%fastapop.strip())
    myfile.close()
def gff2tofu(x,y):
    myfile=open(r'embossORF-%s'%x,'w')#
    for line in open(r'%s'%x):
        if 'CDS\t' in line:
            a=line.split('ID')[0]#注释前半部分
            tranid=(line.split('=')[1]).split('.m')[0]
            geneid='PB.%s'%(tranid.split('.')[1])
            #print '%sgene_id "%s"; transcript_id "%s";\n'%(a,geneid,tranid)
            myfile.write('%sgene_id "%s"; transcript_id "%s";\n'%(a,geneid,tranid))
    for line in open(r'%s'%y):
        myfile.write(line)
    myfile.close()
def FGRAMchange(x,y):#将FGRAM_变成FGRAM.
    myfile=open(r'%s'%y,'w')
    for line1 in open(r'%s'%x):
        if 'FGRAMPH1_' in line1:
            #print line1.replace('FGRAMPH1_','FGRAMPH1.')
            myfile.write(line1.replace('FGRAMPH1_','FGRAMPH1.'))
        if 'FGSG_' in line1:
            myfile.write(line1.replace('FGSG_','FGSG.'))
        else:
            myfile.write(line1)
    myfile.close()
#FGRAMchange('join3-EMBOSS-replaceID.cds','join4-EMBOSS-replaceID.cds')
def RechangeFG(x,y):
    myfile=open(r'%s'%y,'w')
    for line1 in open(r'%s'%x):
        if 'FGRAMPH1.' in line1:
            #print line1.replace('FGRAMPH1_','FGRAMPH1.')
            myfile.write(line1.replace('FGRAMPH1.','FGRAMPH1_'))
        if 'FGSG.' in line1:
            myfile.write(line1.replace('FGSG.','FGSG_'))
        else:
            myfile.write(line1)
    myfile.close()
#FGRAMchange('join3-EMBOSS-replaceID.cds','join4-EMBOSS-replaceID.cds')
def EMBOSSorf(x):
    #1.进行orf预测
    os.system('getorf -table 0 -minsize 300 -find 3 -[no]methionine -reverse false -circular N -flanking 0 -sequence %s -outseq EMBOSS-%s.cds'%(x,x.split('-')[1]))
    #tranid=['>PB.7096.483']
    #2.计算选取最长orf
    joinfasta('EMBOSS-%s.cds'%(x.split('-')[1]),'join-EMBOSS-%s.cds'%(x.split('-')[1]))
    join2fasta('join-EMBOSS-%s.cds'%(x.split('-')[1]),'join2-EMBOSS-%s.cds'%(x.split('-')[1]))
    derev('join2-EMBOSS-%s.cds'%(x.split('-')[1]),'join3-EMBOSS-%s.cds'%((x.split('-')[1])))
    tranid=[]
    for line in open(r'join3-EMBOSS-%s.cds'%x.split('-')[1]):
        if '>' in line:
            tranid.append(line.split('_')[0])
    tranid=list(set(tranid))
    #tranid=['>PB.2949.2']
    longORF=[]
    while len(tranid)>0:
        print len(tranid)
        tranidpop=tranid.pop();tranlen=[]
        for line in open(r'join3-EMBOSS-%s.cds'%x.split('-')[1]):
            if '%s_'%tranidpop in line:
                tranlen.append(int(((line.split('[')[1]).split(']')[0]).split(' - ')[1])-int(((line.split('[')[1]).split(']')[0]).split(' - ')[0]))
        longORF.append('%s\t%s'%(tranidpop,max(tranlen)))
    #3.进行最长ORF的cds提取   
    title=[]
    fasta=[]
    myfile=open(r'%s-longORF.cds'%(x.split('-')[1]),'w')
    for line in open(r'join3-EMBOSS-%s.cds'%(x.split('-')[1])):
        if '>' in line:
            title.append(line)
        else:
            fasta.append(line)
    #print title
    while len(title) and len(fasta) >0:
        #print len(title)
        titlepop=title.pop();fastapop=fasta.pop()
        #print titlepop
        #print fastapop 
        if '%s\t%s'%(titlepop.split('_')[0],int(((titlepop.split('[')[1]).split(']')[0]).split(' - ')[1])-int(((titlepop.split('[')[1]).split(']')[0]).split(' - ')[0])) in longORF:
            myfile.write('%s\n'%titlepop.split('_')[0])
            myfile.write(fastapop)
    myfile.close()
    os.system('/disk/luping/tools/gmap-2017-11-15/bin/gmap -D /disk/luping/bam/iso-bam/pb-cluster-hqlq-total/two-condition/ref-FG -d PH1 -f 2 -n 1 -t 15 --no-chimeras --min-intronlength 20 --max-intronlength-middle 4000 --max-intronlength-ends 4000 -z sense_force %s-longORF.cds>%s-longORF-EMBOSS.gff'%(x.split('-')[1],x.split('-')[1]))
    myfile=open(r'embossORF-%s.gff'%x.split('.')[0],'w')#
    for line in open(r'%s-longORF-EMBOSS.gff'%x.split('-')[1]):
        if 'CDS\t' in line:
            a=line.split('ID')[0]#注释前半部分
            tranid=(line.split('=')[1]).split('.m')[0]
            #geneid='PB.%s'%(tranid.split('.')[1])
            #print '%sgene_id "%s"; transcript_id "%s";\n'%(a,geneid,tranid)
            #myfile.write('%sgene_id "%s"; transcript_id "%s";\n'%(a,geneid,tranid))
            myfile.write('%stranscript_id "%s";\n'%(a,tranid))
    myfile.close()
##首先进行ORF预测
def EmbossORF(x):
    #1 提取所有序列碱基进行ORF预测
    os.system('gffread -w gene-%s.fasta -g ph1.fasta %s'%(x,x))
    #2 ORF 预测
    EMBOSSorf('gene-%s.fasta'%x)
#EmbossORF('PB-replaceID-ReFG-StringTieNew-StingTie1To1-SuppShort-add-Del-classcodeS-delEXP-delLess200bp-C-C-FilterPT-ISO-fusion.collapsed.gff')
#EmbossORF('PB-replaceID-ReFG-StringTieNew-StingTie1To1-SuppShort-add-Del-classcodeS-delEXP-delLess200bp-C-C-FilterPT-ISO-nofusion.collapsed.gff')
#Hy-fusion
EmbossORF('delclassS-C-C-FilterPT-Hy-fusion.collapsed.gff')
#Hy-nofusion
EmbossORF('delclassS-C-C-FilterPT-Hy-nofusion.collapsed.gff')
#Sex-fusion
EmbossORF('delclassS-C-C-FilterPT-Sex-fusion.collapsed.gff')
#Sex-nofusion
EmbossORF('delclassS-C-C-FilterPT-Sex-nofusion.collapsed.gff')
##进行FG的CDS补充
def FGCDS(x,y,z):
    #将PB补充的FG和单独补充的FG转录本号分开
    myfile1=open(r'FGpart-%s'%x,'w')
    PB=[];FG=[]
    for line1 in open(r'%s'%x):
        if 'transcript\t' in line1:
            if 'PB' in line1:
                PB.append( line1.split('"')[3])
            else:
                FG.append( line1.split('"')[3])
    #print len(PB),len(FG)
    for i in PB:#添加PB.FG
        PBFG=i.split('.')[-1]
        for line2 in open(r'%s'%y):
            if 'exon\t' in line2:
                if '=%s.m'%PBFG in line2:
                    CHR=line2.split('\t')[0]
                    source=line2.split('\t')[1]
                    start=line2.split('\t')[3]
                    stop=line2.split('\t')[4]
                    lian1=line2.split('\t')[5]
                    lian2=line2.split('\t')[6]
                    lian3=line2.split('\t')[7]
                    tranid=(line2.split('=')[1]).split('.m')[0]
                    #print '%s\t%s\tCDS\t%s\t%s\t%s\t%s\t%s\ttranscript_id "%s";\n'%(CHR,source,start,stop,lian1,lian2,lian3,i)
                    myfile1.write('%s\t%s\tCDS\t%s\t%s\t%s\t%s\t%s\ttranscript_id "%s";\n'%(CHR,source,start,stop,lian1,lian2,lian3,i))
    for j in FG:#添加FGTAM
        for line3 in open(r'%s'%y):#添加FGRAM
            if 'exon\t' in line3:
                if '=%s.m'%j in line3:
                    CHR=line3.split('\t')[0]
                    source=line3.split('\t')[1]
                    start=line3.split('\t')[3]
                    stop=line3.split('\t')[4]
                    lian1=line3.split('\t')[5]
                    lian2=line3.split('\t')[6]
                    lian3=line3.split('\t')[7]
                    tranid=(line3.split('=')[1]).split('.m')[0]
                    #print '%s\t%s\tCDS\t%s\t%s\t%s\t%s\t%s\ttranscript_id "%s";\n'%(CHR,source,start,stop,lian1,lian2,lian3,j)
                    myfile1.write('%s\t%s\tCDS\t%s\t%s\t%s\t%s\t%s\ttranscript_id "%s";\n'%(CHR,source,start,stop,lian1,lian2,lian3,j))

        for line4 in open(r'%s'%z):#添加FGSG
            if 'exon\t' in line4:
                if '=%sT0'%j in line4:
                    CHR=line4.split('\t')[0]
                    source=line4.split('\t')[1]
                    start=line4.split('\t')[3]
                    stop=line4.split('\t')[4]
                    lian1=line4.split('\t')[5]
                    lian2=line4.split('\t')[6]
                    lian3=line4.split('\t')[7]
                    #print '%s\t%s\tCDS\t%s\t%s\t%s\t%s\t%s\ttranscript_id "%s";\n'%(CHR,source,start,stop,lian1,lian2,lian3,j)
                    myfile1.write('%s\t%s\tCDS\t%s\t%s\t%s\t%s\t%s\ttranscript_id "%s";\n'%(CHR,source,start,stop,lian1,lian2,lian3,j))
    myfile1.close()
#FGCDS('FG-replaceID-ReFG-StringTieNew-StingTie1To1-SuppShort-add-Del-classcodeS-delEXP-delLess200bp-C-C-FilterPT-ISO-fusion.collapsed.gff','gff-FGRAMPH1.gff3','gff-FGSG.gff3')
#FGCDS('FG-replaceID-ReFG-StringTieNew-StingTie1To1-SuppShort-add-Del-classcodeS-delEXP-delLess200bp-C-C-FilterPT-ISO-nofusion.collapsed.gff','gff-FGRAMPH1.gff3','gff-FGSG.gff3')

#接下来，进行gtf二代表达量计算
'''
stringtie /disk/luping/bam/iso-bam/pb-cluster-hqlq-total/two-condition/ISOgene-part/exp-RemainFG/CMC.bam -G PB-FG-ORF-fusion.gff -p 15 -o CMC-fusion.gff
stringtie /disk/luping/bam/iso-bam/pb-cluster-hqlq-total/two-condition/ISOgene-part/exp-RemainFG/YEPD.bam -G PB-FG-ORF-fusion.gff -p 15 -o YEPD-fusion.gff
stringtie /disk/luping/bam/iso-bam/pb-cluster-hqlq-total/two-condition/ISOgene-part/exp-RemainFG/PDA.bam -G PB-FG-ORF-fusion.gff -p 15 -o PDA-fusion.gff
stringtie /disk/luping/bam/iso-bam/pb-cluster-hqlq-total/two-condition/ISOgene-part/exp-RemainFG/CA.bam -G PB-FG-ORF-fusion.gff -p 15 -o CA-fusion.gff
stringtie /disk/luping/bam/iso-bam/pb-cluster-hqlq-total/two-condition/ISOgene-part/exp-RemainFG/TBI.bam -G PB-FG-ORF-fusion.gff -p 15 -o TBI-fusion.gff
stringtie /disk/luping/bam/iso-bam/pb-cluster-hqlq-total/two-condition/ISOgene-part/exp-RemainFG/P6-1.bam -G PB-FG-ORF-fusion.gff -p 15 -o P6-1-fusion.gff
stringtie /disk/luping/bam/iso-bam/pb-cluster-hqlq-total/two-condition/ISOgene-part/exp-RemainFG/P6-2.bam -G PB-FG-ORF-fusion.gff -p 15 -o P6-2-fusion.gff

stringtie /disk/luping/bam/iso-bam/pb-cluster-hqlq-total/two-condition/ISOgene-part/exp-RemainFG/CMC.bam -G PB-FG-ORF-nofusion.gff -p 15 -o CMC-nofusion.gff
stringtie /disk/luping/bam/iso-bam/pb-cluster-hqlq-total/two-condition/ISOgene-part/exp-RemainFG/YEPD.bam -G PB-FG-ORF-nofusion.gff -p 15 -o YEPD-nofusion.gff
stringtie /disk/luping/bam/iso-bam/pb-cluster-hqlq-total/two-condition/ISOgene-part/exp-RemainFG/PDA.bam -G PB-FG-ORF-nofusion.gff -p 15 -o PDA-nofusion.gff
stringtie /disk/luping/bam/iso-bam/pb-cluster-hqlq-total/two-condition/ISOgene-part/exp-RemainFG/CA.bam -G PB-FG-ORF-nofusion.gff -p 15 -o CA-nofusion.gff
stringtie /disk/luping/bam/iso-bam/pb-cluster-hqlq-total/two-condition/ISOgene-part/exp-RemainFG/TBI.bam -G PB-FG-ORF-nofusion.gff -p 15 -o TBI-nofusion.gff
stringtie /disk/luping/bam/iso-bam/pb-cluster-hqlq-total/two-condition/ISOgene-part/exp-RemainFG/P6-1.bam -G PB-FG-ORF-nofusion.gff -p 15 -o P6-1-nofusion.gff
stringtie /disk/luping/bam/iso-bam/pb-cluster-hqlq-total/two-condition/ISOgene-part/exp-RemainFG/P6-2.bam -G PB-FG-ORF-nofusion.gff -p 15 -o P6-2-nofusion.gff
'''
