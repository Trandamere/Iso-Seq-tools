#coding:utf-8
#2020.8.1
#part1,先过滤一波不符合标准的
def NMDfind3(iso):
    #统计有2个CDS转录本的基因
    #part1#找出有两个CDS的转录本
    myfile1=open('NMD-%s'%iso,'w')
    myfile1.write('RefTranscript\tNMD\n')
    Allgene=[]
    for line1 in open(iso):
        if 'transcript\t' in line1:
            Allgene.append( line1.split('gene_id "')[1].split('"')[0])
    Allgene=list(set(Allgene))
    Allgene.sort()
    #统计有两个以上包含CDS转录本的
    for i in Allgene:
        ThaveCDS=[]
        for line2 in open(iso):
            if '%s"'%i in line2:
                if 'PacBio\t' in line2:#三代的外显子
                    if 'CDS\t' in line2:
                        ThaveCDS.append(line2.split('transcript_id "')[1].split('"')[0])
        ThaveCDS=list(set(ThaveCDS))
        if len(ThaveCDS)>1:
            #print i,ThaveCDS
            Texon=[];TCDS=[];Texonlist=[];TCDSlist=[];Strand=[]
            #myfile1.write('%s\n'%i)
            for j in ThaveCDS:
                jexon=[];jCDS=[]
                for line3 in open(iso):
                    if '%s"'%j in line3:
                        if 'exon\t' in line3:
                            jexon.append(int(line3.split('\t')[3]))
                            jexon.append(int(line3.split('\t')[4]))
                            Strand.append(line3.split('\t')[6])
                        if 'CDS\t' in line3:
                            jCDS.append(int(line3.split('\t')[3]))
                            jCDS.append(int(line3.split('\t')[4]))
                jexon.sort()
                Texonlist.append(jexon)
                jCDS.sort()
                TCDSlist.append(jCDS)
                jRange='%s-%s'%(min(jCDS),max(jCDS))
                TCDS.append(jRange)
                #print i,j,jexon,jCDS
                Texon.append(len(jexon))
            TCDS2=list(set(TCDS))
            if Texon.count(2)!=len(Texon):#过滤结果中都是有CDS的单外显子
                if len(TCDS2)!=1:#cds起始终止一样的话，没有PTC
                    #print i,ThaveCDS,Texon,TCDS,TCDS2
                    #首先比较两个转录本CDS起始终止，相同则跳过
                    #然后二者末端外显子是否都有CDS
                    #排除转录起始是否包含CDS
                    #print i,ThaveCDS,Texonlist,TCDSlist
                    Tall=zip(ThaveCDS,Texonlist,TCDSlist)
                    #print Tall
                    Already=[]
                    for J in Tall:
                        for K in Tall:
                            if J[0]!=K[0]:#不和自己比
                                #print J[0],K[0]
                                if Strand[0]=='+':
                                    if J[2][-1]!=K[2][-1]:#排除末端外显子相同
                                        #print J,K
                                        #排除起始终止影响CDS
                                        if J[1][0]<=K[2][0] and J[1][-1]>=K[2][-1] and K[1][0]<=J[2][0] and K[1][-1]>=J[2][-1]:#起始终止位点超过CDS起始终止位点
                                            #print J[0],K[0],'+'
                                            #进行PTC判断，首先要满足有内含子将，即保证有内含子判断
                                            if len(J[1])>2:
                                                LastExonSite=J[1][-3]#末端外显子JS#此处包括两类ixng，自身SJ的CDS在上游50bp，而比较的CDS位点长于SJ，和相反的
                                                #自身CDS长于SJ而比较的CDS位点在SJ上游50bp以上位置
                                                if LastExonSite-J[2][-1]>=50:#自身末端外显子链接位点距离自身CDS超过50bp#正链CDS值小
                                                   if LastExonSite<=K[2][-1]:#比较的转录本，CDS终止位点在SJ位点之后#正链，比较CDS靠后，数值大
                                                       #print J[0],K[0],'+'
                                                       if '%s-%s'%(K[0],J[0]) not in Already:
                                                           print K[0],'exon',J[0],'NMD','+'
                                                           myfile1.write('%s\t%s\n'%(K[0],J[0]))
                                                           Already.append('%s-%s'%(K[0],J[0]))
                                                if LastExonSite<=J[2][-1]:#z自身CDS终止位点要长于SJ，#但比较转录本，CDS位点在SJ位点上游50bp#CDS大，靠后有
                                                   if LastExonSite-K[2][-1]>=50:#C符合缩短，CDS靠前，要小
                                                       if '%s-%s'%(J[0],K[0]) not in Already:
                                                           print J[0],'exon',K[0],'NMD','+'
                                                           myfile1.write('%s\t%s\n'%(J[0],K[0]))
                                                           Already.append('%s-%s'%(J[0],K[0]))
                                if Strand[0]=='-':
                                    if J[2][0]!=K[2][0]:#排除末端外显子相同
                                        #print J[0],K[0],'-'
                                        if J[1][0]<=K[2][0] and J[1][-1]>=K[2][-1] and K[1][0]<=J[2][0] and K[1][-1]>=J[2][-1]:#起始终止位点超过CDS起始终止位点
                                            if len(J[1])>2:
                                                LastExonSite=J[1][2]#负链，的末端外显子JS
                                                #print 'LastExonSite',LastExonSite
                                                if J[2][0]-LastExonSite>=50:#自身末端外显子链接位点距离自身CDS超过50bp#负链，CDS大
                                                    #print J[0],K[0],'more 50'
                                                    if LastExonSite>=K[2][0]:#比较的转录本，CDS终止位点在SJ位点之后#CDS长，则CDS小
                                                        if '%s-%s'%(K[0],J[0]) not in Already:
                                                            print K[0],'exon',J[0],'NMD','-'
                                                            myfile1.write('%s\t%s\n'%(K[0],J[0]))
                                                            Already.append('%s-%s'%(K[0],J[0]))
                                                if J[2][0]<=LastExonSite:#自身CDS终止位点要长于SJ，#但比较转录本，CDS位点在SJ位点上游50bp#在前方，前方小#CDS小
                                                    if K[2][0]-LastExonSite>=50:#CDS短，CDS大
                                                        if '%s-%s'%(J[0],K[0]) not in Already:
                                                            myfile1.write('%s\t%s\n'%(J[0],K[0]))
                                                            print J[0],'exon',K[0],'NMD','-'
                                                            Already.append('%s-%s'%(J[0],K[0]))

    myfile1.close()

#NMDfind3('PH1-fusion-version3.3-iso.gtf')

#NMDfind3('Hy-fusion-FixJS.gtf')
#NMDfind3('Sex-fusion-FixJS.gtf')
#NMDfind3('iso-PH1-fusion-version4.1.gtf')#修正NMD筛选

#进一步筛选去除，Ref不是参考的转录本
def FilNoRef(x,y):
    myfile1=open('MaxT-ref-%s'%y,'w')
    myfile2=open('MaxT-noref-%s'%y,'w')
    #T=[];F=[]
    Ref=[]
    for line1 in open(x):
        Ref.append( line1.strip().split('\t')[13])
    for line2 in open(y):
        if line2.split('\t')[0] in Ref:
            myfile1.write(line2)
            #T.append(line2.split())
        else:
            myfile2.write(line2)
    myfile1.close()
    myfile2.close()
#FilNoRef('Gene-T-for-Allgene.txt','NMD-PH1-fusion-version3.3-iso.gtf')
FilNoRef('Gene-T-for-Allgene.txt','NMD-iso-PH1-fusion-version4.1.gtf')
