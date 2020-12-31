#coding:utf-8
#2020.12.31

#FLNC的PAS位点上下哟
def delNA4(PASfile,gff):
    myfile=open('PAS-updown201bp.gtf','w')
    ID=[]
    for line1 in open(PASfile):
        ID.append( line1.split('\t')[0])
    ID=list(set(ID));ID.sort()
    #print len(ID)
    for i in ID:
        print i
        PAS=[];Strand=[];CHR=[]
        for line2 in open(PASfile):
            if '%s\t'%i in line2:
                PAS.append(int(line2.strip().split('-')[1]))
        for line3 in open(gff):
            if '%s"'%i in line3:
                Strand.append(line3.split('\t')[6])
                CHR.append(line3.split('\t')[0])
        #print i,PAS,Strand[0]
        for k in PAS:
            #print '%s\tPacBio\texon\t%s\t%s\t.\t%s\t.\tgene_id "%s"; transcript_id "%s-%s";\n'%(CHR[0],k-100,k+100,Strand[0],i,i,k)
            myfile.write('%s\tPacBio\texon\t%s\t%s\t.\t%s\t.\tgene_id "%s"; transcript_id "%s-%s";\n'%(CHR[0],k-100,k+100,Strand[0],i,i,k))
    myfile.close()

#delNA4('All-PAS-number-count.txt','PH1-fusion-version4.2.gtf')
#gffread -w PAS-updown201bp-all.fasta -g ph1.fasta PAS-updown201bp.gtf


#fasta序列展开
def joinfasta(x,y):
    myfile=open(r'%s'%y,'w')
    for line in open(r'%s'%x):
        if '>' in line:
            #print line
            myfile.write('\n%s'%line)
        else:
            #print line.strip()
            myfile.write(line.strip())
    myfile.close()    
#joinfasta('PAS-updown201bp-all.fasta','join-PAS-updown201bp-all.fasta')

##将所有碱基分开成
def Nfrequency(fasta):
    myfile=open('%s-GTAG-frequency.txt'%fasta,'w')
    for line1 in open(fasta):
        if '>' not in line1:
            FG='%s'%(list(line1.strip()))
            #print FG.replace("['","").replace("', '","\t").replace("']","\n")
            myfile.write(FG.replace("['","").replace("', '","\t").replace("']","\n"))
    myfile.close()
#Nfrequency('join-PAS-updown201bp-all.fasta')


#统计201个碱基上的CTAG百分比
def CTAGf(x):
    myfile=open('Count-%s'%x,'w')
    #n=0;m=-50#200bp替换处
    n=0;m=-100
    myfile.write('Position\tDNA\tFrequency\n')
    #while n<=100:#200bp替换处
    while n<=200:
        print n
        FGA=[];FGT=[];FGC=[];FGG=[]
        for line1 in open(x):
            #print m,line1.strip().split('\t')[n]
            bp=line1.strip().split('\t')[n]
            if bp=='A':
                FGA.append(bp)
            if bp=='T':
                FGT.append(bp)
            if bp=='C':
                FGC.append(bp)
            if bp=='G':
                FGG.append(bp)
        Total=len(FGA)+len(FGT)+len(FGC)+len(FGG)
        #print m,'A',float(len(FGA))/float(Total)
        #print m,'T',float(len(FGT))/float(Total)
        #print m,'C',float(len(FGC))/float(Total)
        #print m,'G',float(len(FGG))/float(Total)
        myfile.write('%s\t%s\t%s\n'%(m,'A',float(len(FGA))/float(Total)))
        myfile.write('%s\t%s\t%s\n'%(m,'T',float(len(FGT))/float(Total)))
        myfile.write('%s\t%s\t%s\n'%(m,'C',float(len(FGC))/float(Total)))
        myfile.write('%s\t%s\t%s\n'%(m,'G',float(len(FGG))/float(Total)))
        n+=1
        m+=1
    myfile.close()
CTAGf('join-PAS-updown201bp-all.fasta-GTAG-frequency.txt')
