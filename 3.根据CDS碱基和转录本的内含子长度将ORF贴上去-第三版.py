#coding:utf-8
import numpy
def FGid(cds,gff):
    myfile=open('%s-CDS.gtf'%cds,'w')
    T=[]#FG1G00820T4
    for line1 in open(gff):
        if 'transcript\t' in line1:
            T.append( line1.split('transcript_id "')[1].split('"')[0])
    T.sort()
    #T=['FG1G00790T1']#有问题的反向时期ATG#FG1G02180T3#完成
    #T=['FG1G00910T1']#3UTR有内含子#完成
    #T=['FG1G07430T0','FG1G07430T4','FG1G07430T5']#完成
    #T=['FG1G01100T0']#起始位点ATG#完成
    #T=['FG1G00420T3']#假APA位点计算#完成
    #T=['FG1G00650T8']#反向APA假位点#完成
    #T=['FG1G10380T1']
    #T=['FG2G05010T6']#只有一个外显子且在最末端
    #T=['FG1G36520T2']#起始位点在第二个内含子#ATG在第二个内含子处
    #T=['FG1G16230T2']#末端外显子补3的差别判断
    
    for i in T:
        print i
        CDSstart=[];CDSL=[]
        for line2 in open(cds):
            if '>%s_'%i in line2:
                #print i,line2#获取cds起始终止区域
                CDSstart.append(int(line2.split(' - ')[0].split('[')[1]))
                CDSL.append(abs(int(line2.split(' - ')[0].split('[')[1])-int(line2.split(' - ')[1].split(']')[0])))
        #print CDSstart,CDSL#接下来获得转录本内含子总长度和链方向
        exon=[];strand=[];exon2=[];exon3=[];Source=[];CHR=[];Gene=[]
        for line3 in open(gff):
            if 'exon\t' in line3:
                if '%s"'%i in line3:
                    exon.append(int(line3.split('\t')[3]))
                    exon.append(int(line3.split('\t')[4]))#exon2.append(int(line3.split('\t')[3]))
                    exon2.append(int(line3.split('\t')[4]))
                    exon3.append(int(line3.split('\t')[3]))#exon3.append(int(line3.split('\t')[4]))
                    strand.append(line3.split('\t')[6])
                    Source.append(line3.split('\t')[1])
                    CHR.append(line3.split('\t')[0])
                    Gene.append(line3.split('gene_id ')[1])
        exon.sort();exon2.sort();exon3.sort()
        del exon2[-1];del exon3[0]
        #del exon2[0];del exon3[-1]
        #print exon2,exon3,strand
        exon2=numpy.array(exon2)
        exon3=numpy.array(exon3)
        Intron=exon3-exon2
        IntronL=sum(Intron)
        #print Intron
        #print 'min(exon)',min(exon),'Intron',Intron,'IntronL',IntronL,'CDSstart',CDSstart,'IntronL',IntronL,'CDSL',CDSL[0]
        if len(Intron)>0:#有内含子#有内含子在起始位点前方的还要计算起始位点加上内含子位点#将内含子分为起始位点前和起始位点后
            if len(CDSstart)>0:#需要考虑ATG起始位点的问题。区分假APA位点计算和，TAG挎在内含子区域区分的问题，来决定删除哪个位点。FG1G07430
                if strand[0]=='+':
                    #找出内含子的插入位置，比较，分类#第一个内含子位置为第二个剪切位点
                    BeforeIntron=[];AfterIntron=[];UTR3Intron=[]#分成起始位点前后的内含子及长度，起始位点计算不同#此外还要计算3UTR的内含子长度，不在计算CDS的范围
                    num=0
                    CDSexon=[]
                    HaveATGintron=[]
                    HaveTAAintron=[]#有终止位点挎在内含子上
                    for N in exon2:
                        #print N
                        if num==0:
                            #print N-min(exon),'CDSstart',CDSstart
                            #print 'Intron index',Intron[num]
                            #print 'HAAH0',(int(N)-min(exon))+1,CDSstart[0]
                            if (int(N)-min(exon))+1<CDSstart[0]:
                                #print int(N)-min(exon)+1,CDSstart[0]
                                BeforeIntron.append(Intron[num])
                                num+=1
                            elif (int(N)-min(exon)+1)==CDSstart[0]:##ATG的A正好在内含子的边界上
                                AfterIntron.append(Intron[num])
                                HaveATGintron.append(min(exon)+CDSstart[0]-1+sum(BeforeIntron)-len(BeforeIntron))
                                num+=1
                            elif (int(N)-min(exon)+1)==(CDSstart[0]+CDSL[0]):
                                AfterIntron.append(Intron[num])
                                HaveTAAintron.append((min(exon)+CDSstart[0]-1+sum(BeforeIntron)-len(BeforeIntron))+(CDSL[0]+sum(AfterIntron)-len(AfterIntron)+3))
                                num+=1
                            elif CDSstart[0]<(int(N)-min(exon)+1)<(CDSstart[0]+CDSL[0]):
                                AfterIntron.append(Intron[num])
                                num+=1
                            else:
                                UTR3Intron.append(Intron[num])
                                num+=1
                        else:
                            #print 'before,intron',sum(Intron[0:num])
                            #print N-min(exon)-sum(Intron[0:num]),'CDSstart',CDSstart
                            #print 'Intron index',Intron[num]
                            if (int(N)-min(exon)+1-sum(BeforeIntron)-sum(AfterIntron)+len(BeforeIntron)+len(AfterIntron))<CDSstart[0]:
                                BeforeIntron.append(Intron[num])
                                num+=1
                            elif (int(N)-min(exon)+1-sum(BeforeIntron)-sum(AfterIntron)+len(BeforeIntron)+len(AfterIntron))==CDSstart[0]:##ATG的A正好在内含子的边界上
                                AfterIntron.append(Intron[num])
                                #HaveATGintron.append(min(exon)+CDSstart[0]-1+sum(BeforeIntron)-len(BeforeIntron))
                                num+=1
                            elif (int(N)-min(exon)+1-sum(BeforeIntron)-sum(AfterIntron)+len(BeforeIntron)+len(AfterIntron))==(CDSstart[0]+CDSL[0]):
                                AfterIntron.append(Intron[num])
                                #HaveTAAintron.append(CDSstartSite+CDSL[0]+sum(AfterIntron)-len(AfterIntron)+3)
                                num+=1
                            elif CDSstart[0]<(int(N)-min(exon)+1-sum(BeforeIntron)-sum(AfterIntron)+len(BeforeIntron)+len(AfterIntron))<(CDSstart[0]+CDSL[0]):
                                AfterIntron.append(Intron[num])
                                num+=1
                            else:
                                UTR3Intron.append(Intron[num])
                                num+=1
                    print 'CDS exon',CDSexon
                    print 'Before,After,UTR3',BeforeIntron,AfterIntron,UTR3Intron
                    print 'HaveATGintron','HaveTAAintron',HaveATGintron,HaveTAAintron
                    CDSstartSite=min(exon)+CDSstart[0]-1+sum(BeforeIntron)-len(BeforeIntron)
                    CDSstopSite=CDSstartSite+CDSL[0]+sum(AfterIntron)-len(AfterIntron)+3#减去内含子带来的长度影响，再加上缺失的终止密码子+1+2
                    print 'CDS start stop',CDSstartSite,CDSstopSite,'!'*50#生成CDS打印的exon区域#直接读取转录本exon区域，生成
                    CDSexon.append(CDSstartSite);CDSexon.append(CDSstopSite)
                    if abs(max(CDSexon)-max(exon))<=3:
                        print 'abs(max(CDSexon)-max(exon))',max(CDSexon),max(exon)
                        del CDSexon[1];CDSexon.append(max(exon))#这里去除了末端外显子与CDS查3的结果
                    print 'CDS exon',CDSexon
                    for j in exon:
                        if CDSstartSite<=j<=CDSstopSite+3:#这里要加上差3的结果
                            print 'j',j
                            CDSexon.append(j)                  
                    print 'CDS exon',CDSexon
                    CDSexon=list(set(CDSexon)) 
                    if len(HaveATGintron)>0:
                        CDSexon.append(HaveATGintron[0])
                    if len(HaveTAAintron)>0:
                        CDSexon.append(HaveTAAintron[0])
                    CDSexon.sort()
                    if len(CDSexon)==1:###只有一个位点的时候添加
                        CDSexon.append(CDSstopSite)
                    CDSexon.sort()
                    print 'test CDS exon',CDSexon
                    if (len(CDSexon)%2)==0:#CDS的位点是双数，说明终止位点在预测范围内，不需要矫正
                        pass
                    else:
                        del CDSexon[-1]#如果为奇数，说明终止位点不在预测范围内，只用原来的位点
                    n=0
                    print 'CDS exon',CDSexon
                    while n<len(CDSexon):
                        #CDS数量是偶数#有终止位点
                        print '%s\t%s\tCDS\t%s\t%s\t.\t%s\t.\tgene_id %s'%(CHR[0],Source[0],CDSexon[n],CDSexon[n+1],strand[0],Gene[0])
                        myfile.write('%s\t%s\tCDS\t%s\t%s\t.\t%s\t.\tgene_id %s'%(CHR[0],Source[0],CDSexon[n],CDSexon[n+1],strand[0],Gene[0]))
                        n+=2
                if strand[0]=='-':
                    Intron2=Intron[::-1]#内含子位置反转
                    #print 'Intron2',Intron2
                    exon4=exon3[::-1]
                    BeforeIntron=[];AfterIntron=[];UTR3Intron=[]#分成起始位点前后的内含子及长度，起始位点计算不同
                    num=0
                    CDSexon=[]
                    HaveATGintron=[]
                    HaveTAAintron=[]#有终止位点挎在内含子上
                    for N in exon4:
                        if num==0:
                            print 'HAAH0',max(exon)-int(N)-sum(BeforeIntron)+len(BeforeIntron)-CDSstart[0]
                            if max(exon)-int(N)+1<CDSstart[0]:#少了加1的计算
                                #print max(exon),int(N),CDSstart
                                BeforeIntron.append(Intron2[num])
                                num+=1
                            elif max(exon)-int(N)+1==CDSstart[0]:##起始位点挎在第一个内含子
                                AfterIntron.append(Intron2[num])
                                #HaveATGintron.append(1)
                                HaveATGintron.append(max(exon)-CDSstart[0]+1-sum(BeforeIntron)+len(BeforeIntron))
                                print 'HaveATG1',max(exon)-CDSstart[0]+1-sum(BeforeIntron)+len(BeforeIntron)
                                num+=1
                            elif max(exon)-int(N)+1==(CDSstart[0]+CDSL[0]):
                                AfterIntron.append(Intron2[num])
                                #HaveTAAintron.append((max(exon)-CDSstart[0]+1-sum(BeforeIntron)+len(BeforeIntron))-(CDSstartSite-CDSL[0]-sum(AfterIntron)+len(AfterIntron)-3))
                                #HaveTAAintron.append((max(exon)-CDSstart[0]+1-sum(BeforeIntron)+len(BeforeIntron))-(CDSL[0]-sum(AfterIntron)+len(AfterIntron)-3))
                                num+=1
                            elif CDSstart[0]<max(exon)-int(N)+1<(CDSstart[0]+CDSL[0]):
                                AfterIntron.append(Intron2[num])
                                num+=1
                            else:
                                UTR3Intron.append(Intron2[num])
                                num+=1
                        else:
                            print 'HAAH>0',max(exon)-int(N)-sum(BeforeIntron)-sum(AfterIntron)+len(BeforeIntron)+len(AfterIntron)-CDSstart[0]
                            if max(exon)-int(N)+1-sum(BeforeIntron)-sum(AfterIntron)+len(BeforeIntron)+len(AfterIntron)<CDSstart[0]:
                                BeforeIntron.append(Intron2[num])
                                num+=1
                            elif max(exon)-int(N)+1-sum(BeforeIntron)-sum(AfterIntron)+len(BeforeIntron)+len(AfterIntron)==CDSstart[0]:##起始位点挎在第一个内含子
                                AfterIntron.append(Intron2[num])
                                #HaveATGintron.append(max(exon)-CDSstart[0]+1-sum(BeforeIntron)+len(BeforeIntron))
                                #HaveATGintron.append((max(exon)-CDSstart[0]+1-sum(BeforeIntron)+len(BeforeIntron))-CDSL[0]-sum(AfterIntron)+len(AfterIntron)-3)
                                num+=1
                            elif max(exon)-int(N)+1-sum(BeforeIntron)-sum(AfterIntron)+len(BeforeIntron)+len(AfterIntron)==(CDSstart[0]+CDSL[0]):
                                AfterIntron.append(Intron2[num])
                                #HaveTAAintron.append((max(exon)-CDSstart[0]+1-sum(BeforeIntron)+len(BeforeIntron))-(CDSstartSite-CDSL[0]-sum(AfterIntron)+len(AfterIntron)-3))
                                #HaveTAAintron.append(max(exon)-CDSstart[0]+1-sum(BeforeIntron)+len(BeforeIntron))
                                num+=1
                            elif CDSstart[0]<max(exon)-int(N)+1-sum(BeforeIntron)-sum(AfterIntron)+len(BeforeIntron)+len(AfterIntron)<(CDSstart[0]+CDSL[0]):
                                AfterIntron.append(Intron2[num])
                                num+=1
                            else:
                                UTR3Intron.append(Intron2[num])
                                num+=1
                    print  'Before,After,UTR3',BeforeIntron,AfterIntron,UTR3Intron
                    CDSstartSite=max(exon)-CDSstart[0]+1-sum(BeforeIntron)+len(BeforeIntron)
                    CDSstopSite=CDSstartSite-CDSL[0]-sum(AfterIntron)+len(AfterIntron)-3
                    print 'CDS start stop',CDSstartSite,CDSstopSite,'!'*50
                    CDSexon.append(CDSstartSite);CDSexon.append(CDSstopSite)
                    print 'CDS exon',CDSexon
                    if abs(min(CDSexon)-min(exon))<=3:
                        del CDSexon[1];CDSexon.append(min(exon))
                    print 'min(CDSexon)-min(exon)',min(CDSexon)-min(exon)
                    print 'CDS exon',CDSexon
                    for k in exon:
                        if CDSstopSite-3<=k<=CDSstartSite:
                            print 'k',k
                            CDSexon.append(k)
                    CDSexon=list(set(CDSexon))
                    CDSexon.sort()
                    print 'CDS exon',CDSexon
                    '''
                    if (len(CDSexon)%2)==0:#CDS的位点是双数，说明终止位点在预测范围内，不需要矫正
                        pass
                    else:
                        #del CDSexon[-1]#如果为奇数，说明终止位点不在预测范围内，只用原来的位点
                        del CDSexon[1]#这个删除查看结果
                        #
                    print 'CDS exon',CDSexon
                    '''
                    if len(HaveATGintron)>0:
                        CDSexon.append(HaveATGintron[0])
                    if len(HaveTAAintron)>0:
                        CDSexon.append(HaveTAAintron[0])
                    CDSexon.sort()
                    print 'CDS exon',CDSexon
                    print 'HaveATGintron','HaveTAAintron',HaveATGintron,HaveTAAintron
                    ####下面这步时为了防止除了ATG起始外，还存在起始终止的位点选择
                    if len(CDSexon)==1:###只有一个位点的时候添加
                        CDSexon.append(CDSstopSite)
                    print 'test CDS exon',CDSexon
                    if len(CDSexon)>2:
                        if (len(CDSexon)%2)==0:#CDS的位点是双数，说明终止位点在预测范围内，不需要矫正
                            pass
                        else:
                            del CDSexon[1]#如果为奇数，说明终止位点不在预测范围内，只用原来的位点
                    CDSexon.sort()
                    n=0
                    while n<len(CDSexon):
                        #CDS数量是偶数#有终止位点
                        print '%s\t%s\tCDS\t%s\t%s\t.\t%s\t.\tgene_id %s'%(CHR[0],Source[0],CDSexon[n],CDSexon[n+1],strand[0],Gene[0])
                        myfile.write('%s\t%s\tCDS\t%s\t%s\t.\t%s\t.\tgene_id %s'%(CHR[0],Source[0],CDSexon[n],CDSexon[n+1],strand[0],Gene[0]))
                        n+=2
        else:#没有内含子的直接计算
            #print 'min(exon)',min(exon),'Intron',Intron,'IntronL',IntronL,'CDSstart',CDSstart,'CDSL',CDSL[0]
            if len(CDSstart)>0:
                if strand[0]=='+':
                    CDSstartSite=min(exon)+CDSstart[0]-1
                    CDSstopSite=CDSstartSite+CDSL[0]+3
                    #print 'CDSstartSite',CDSstartSite,'CDSstopSite',CDSstopSite
                    #print '%s\t%s\tCDS\t%s\t%s\t.\t%s\t.\tgene_id %s'%(CHR[0],Source[0],CDSstartSite,CDSstopSite,strand[0],Gene[0])
                    myfile.write('%s\t%s\tCDS\t%s\t%s\t.\t%s\t.\tgene_id %s'%(CHR[0],Source[0],CDSstartSite,CDSstopSite,strand[0],Gene[0]))
                if strand[0]=='-':
                    CDSstartSite=max(exon)-CDSstart[0]+1
                    CDSstopSite=CDSstartSite-CDSL[0]-3
                    #print 'CDSstartSite',CDSstopSite,'CDSstopSite',CDSstartSite
                    #print '%s\t%s\tCDS\t%s\t%s\t.\t%s\t.\tgene_id %s'%(CHR[0],Source[0],CDSstopSite,CDSstartSite,strand[0],Gene[0])
                    myfile.write('%s\t%s\tCDS\t%s\t%s\t.\t%s\t.\tgene_id %s'%(CHR[0],Source[0],CDSstopSite,CDSstartSite,strand[0],Gene[0]))
    myfile.close()
FGid('Fusion-embossORF-gene-delclassS-C-C-FilterPT-Hy-fusion.gff.gff3.gtf-longORF.cds','Fusion-Hy-delMt.gtf')
