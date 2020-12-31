#coding:utf-8
def FGid(PAScountFile,ClassTfile,UTRfile,Intronfile):
    myfile=open('PAS-location-count.txt','w')
    myfile.write('Gene\tClass Transcript\tPASs\tPAS\tPAS Type\tStrand\tCDS range\t5UTR range\t3UTR range\tIntron range\tTranscript range\n')
    ID=[]
    for line1 in open(PAScountFile):
        ID.append(line1.split('\t')[0])
    ID=list(set(ID));ID.sort()
    for i in ID:
        print i
        ClassT=[];UTR3=[];UTR5=[];CDS=[];Intron=[];Exon=[];Strand=[];Trange=[];PAS=[];UTR3range=[];UTR5range=[];CDSrange=[];Intronrange=[]
        for line2 in open(ClassTfile):
            if '%s\t'%i in line2:#print i,line2.strip().split('\t')[13]
                ClassT.append(line2.strip().split('\t')[13])#print i,ClassT
        for line3 in open(UTRfile):#收集链方向，CDS区域，5UTR和3UTR，转录本的位点。
            if '%s;'%ClassT[0] in line3:
                Strand.append(line3.split('\t')[6])
                if 'CDS\t' in line3:
                    CDSrange.append('%s-%s'%(int(line3.split('\t')[3]),int(line3.split('\t')[4])))
                    for I in range(int(line3.split('\t')[3]),int(line3.split('\t')[4])+1):
                        CDS.append(I)
                if 'five_prime_UTR\t' in line3:
                    UTR5range.append('%s-%s'%(int(line3.split('\t')[3]),int(line3.split('\t')[4])))
                    for J in range(int(line3.split('\t')[3]),int(line3.split('\t')[4])+1):
                        UTR5.append(J)
                if 'three_prime_UTR\t' in line3:
                    UTR3range.append('%s-%s'%(int(line3.split('\t')[3]),int(line3.split('\t')[4])))
                    for K in range(int(line3.split('\t')[3]),int(line3.split('\t')[4])+1):
                        UTR3.append(K)
                if 'mRNA\t' in line3:
                    Trange.append(int(line3.split('\t')[3]))
                    Trange.append(int(line3.split('\t')[4]))
        for line4 in open(Intronfile):
            if '%s-'%ClassT[0] in line4:
                Intronrange.append('%s-%s'%(int(line3.split('\t')[3]),int(line3.split('\t')[4])))
                for L in range(int(line4.split('\t')[3]),int(line4.split('\t')[4])+1):
                    Intron.append(L)
        for line5 in open(PAScountFile):
            if '%s\t'%i in line5:
                PAS.append(int(line5.strip().split('-')[1]))
        #此处开始判断位点所处的区域#分为有CDS区域和没有CDS区域的。有CDS区域的可进行CDSUTR，inron判定。
        if len(CDS)>0:#有CDS的，有UTR
            for PASi in PAS:
                #print i,PAS,PASi
                if PASi in UTR3:
                    #print i,ClassT[0],PASi,'3UTR'
                    myfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(i,ClassT[0],PAS,PASi,'3UTR',Strand[0],CDSrange,UTR5range,UTR3range,Intronrange,Trange))
                elif PASi in CDS:
                    #print i,ClassT[0],PASi,'CDS'
                    myfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(i,ClassT[0],PAS,PASi,'CDS',Strand[0],CDSrange,UTR5range,UTR3range,Intronrange,Trange))
                elif PASi in UTR5:
                    #print i,ClassT[0],PASi,'5UTR'
                    myfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(i,ClassT[0],PAS,PASi,'5UTR',Strand[0],CDSrange,UTR5range,UTR3range,Intronrange,Trange))
                elif PASi in Intron:
                    #print i,ClassT[0],PASi,'Intron'
                    myfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(i,ClassT[0],PAS,PASi,'Intron',Strand[0],CDSrange,UTR5range,UTR3range,Intronrange,Trange))
                else:
                    #print i,ClassT[0],PASi,Trange,'What?'*100#超出主要转录本范围，按照转录本链方向判定5UTR和3UTR
                    if Strand[0]=='+':
                        if max(Trange)<PASi:
                            #print i,ClassT[0],PASi,'3UTR'
                            myfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(i,ClassT[0],PAS,PASi,'3UTR',Strand[0],CDSrange,UTR5range,UTR3range,Intronrange,Trange))
                        if PASi<min(Trange):
                            #print i,ClassT[0],PASi,'5UTR'
                            myfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(i,ClassT[0],PAS,PASi,'5UTR',Strand[0],CDSrange,UTR5range,UTR3range,Intronrange,Trange))
                    if Strand[0]=='-':
                        if PASi<min(Trange):
                            #print i,ClassT[0],PASi,'3UTR'
                            myfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(i,ClassT[0],PAS,PASi,'3UTR',Strand[0],CDSrange,UTR5range,UTR3range,Intronrange,Trange))
                        if max(Trange)<PASi:
                            #print i,ClassT[0],PASi,'5UTR'
                            myfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(i,ClassT[0],PAS,PASi,'5UTR',Strand[0],CDSrange,UTR5range,UTR3range,Intronrange,Trange))
        else:#没有CDS的，只有外显子和内含子之分。
            if PASi in Intron:
                #print i,ClassT[0],PASi,'Intron'
                myfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(i,ClassT[0],PAS,PASi,'Intron',Strand[0],CDSrange,UTR5range,UTR3range,Intronrange,Trange))
            else:#没有CDS只有内含子和外显子，不在内含子，就是在外显子上。
                #print i,ClassT[0],PASi,'Exon'
                myfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(i,ClassT[0],PAS,PASi,'Exon',Strand[0],CDSrange,UTR5range,UTR3range,Intronrange,Trange))
    myfile.close()

#FGid('All-PAS-number-count.txt','Gene-T-for-Allgene-add-noCDS.txt','PH1-fusion-version4.4.gff3','Fusion-PH1-fusion-version2.3.gtf.gff3.gtf-tran-intron.gff')
def CountFLNCnum(x,y):
    myfile=open('add-FLNCcount-%s'%x,'w')
    for line1 in open(x):
        if 'FG' not in line1:
            #print '%s\tFLNCcount\n'%line1.strip()
            myfile.write('%s\tFLNCcount\n'%line1.strip())
        else:
            Site= '%s-%s'%(line1.split('\t')[0],line1.split('\t')[3])
            print Site
            Num=[]
            for line2 in open(y):
                if '%s\t'%Site in line2:
                    Num.append(int(line2.strip().split('\t')[2]))
            if len(Num)>0:#有数量统计的为APA部分
                #print '%s\t%s\n'%(line1.strip(),Num[0])
                myfile.write('%s\t%s\n'%(line1.strip(),Num[0]))
            else:#没有数量统计的，为最初统计为单个PAS序列的
                #print '%s\t%s\n'%(line1.strip(),1)
                myfile.write('%s\t%s\n'%(line1.strip(),1))
CountFLNCnum('PAS-location-count.txt','All-APA-PAS-Number-count.txt')
#统计文件附加FLNC的数量统计，以计算非uniq的数量
