#coding:utf-8
#选取PacBio基因的主要表达转录本为结构，生成roar注释
class Solution:
    def removeOuterParentheses(self, S):
        stack = []
        label = []
        string = ''
        for i in range(len(S)):
            if len(stack)==0: label.append(i) # mark the start
            if stack==[]:
                stack.append(S[i])
            else:
                if S[i]=='[': stack.append(S[i])
                if S[i]==']': stack.pop()
            if len(stack)==0: label.append(i+1) # mark the end
        for i in range(len(label)):
            if i%2==0:
                string += S[label[i]:label[i+1]][1:-1]
        return string

def roar2(count,countNum,Ref,gff):#修正，去除PAS位点在主转录本起始位点前方的位点
    myfile=open('roar2.gtf','w')
    for line1 in open(count):
        gene=[];CHR=[];strand=[];APA=[];site=[];transcript=[]
        if 'FG' in line1:
            if int(line1.split('\t')[3])>1: #APA数量大于1
                gene.append(line1.split('\t')[0])
                CHR.append(line1.split('\t')[1])
                strand.append(line1.split('\t')[2])
                APA.append(int(line1.split('\t')[3]))
                poly=line1.split('\t')[4].replace('"',"")
                #print poly
                Location=Solution().removeOuterParentheses(poly)#去掉位点的外层方括号
                #print Location
                Location2= Location.split(']')
                #print Location2
                for i in Location2:
                        #print 'i',i
                        if ', [' in i:#接下来区分有多个位点和单个位点的列表，单个位点直接添加，多个位点的经过文件验证选取支持位点最多的位点
                            PASlist=i.split(', [')
                            PASlist.remove('')
                            #print gene[0],PASlist,len(PASlist),'list1'#选取位点的第一个
                            #print len(PASlist)
                            #site.append(int(i.split(', [')[1].split(',')[0]))
                            if ',' not in PASlist[0]:#没有逗号说明只有一个位点
                                site.append(int(PASlist[0]))
                            else:
                                PASnum=[];PASlocationlist1=[]
                                MpasList= PASlist[0].split(', ')
                                for j in MpasList:
                                    GeneSite= '%s-%s'%(gene[0],j)
                                    for line2 in open(countNum):
                                        if '%s\t'%GeneSite in line2:
                                            #print j,GeneSite,line2,int(line2.strip().split('\t')[2])
                                            PASlocationlist1.append(j)
                                            PASnum.append(int(line2.strip().split('\t')[2]))
                                #print PASlocationlist1,PASnum,max(PASnum),PASlocationlist1[PASnum.index(max(PASnum))]
                                site.append(int(PASlocationlist1[PASnum.index(max(PASnum))]))#选取支持位点数最大的转录本                                         
                        elif '[' in i:
                            PASlist2=i.split('[')
                            PASlist2.remove('')
                            #print gene[0],PASlist2,len(PASlist2),'list2'#print len(PASlist2)
                            if ',' not in PASlist2[0]:#没有逗号说明只有一个位点
                                site.append(int(PASlist2[0]))
                            else:
                                PASnum2=[];PASlocationlist2=[]
                                MpasList2= PASlist2[0].split(', ')
                                for J in MpasList2:
                                    GeneSite2= '%s-%s'%(gene[0],J)
                                    for line3 in open(countNum):
                                        if '%s\t'%GeneSite2 in line3:
                                            #print J,GeneSite2,line3,int(line3.strip().split('\t')[2])
                                            PASlocationlist2.append(J)
                                            PASnum2.append(int(line3.strip().split('\t')[2]))
                                site.append(int(PASlocationlist2[PASnum2.index(max(PASnum2))]))#选取支持位点数最大的转录本    
        #print gene,CHR,strand,APA,transcript,site
        if len(gene)>0:#找出基因主要表达转录本
            print gene
            if strand[0]=='+':
                #MinPAS=min(site);MaxPAS=max(site)
                #print gene,MinPAS,MaxPAS
                for I in site:
                    #print '%s\tPacBio\tapa\t%s\t%s\t.\t+\t.\tapa "%s.%s_%s"\n'%(CHR[0],I,I,site.index(I),I,gene[0])
                    myfile.write('%s\tPacBio\tapa\t%s\t%s\t.\t+\t.\tapa "%s.%s_%s"\n'%(CHR[0],I,I,site.index(I),I,gene[0]))
                #编写基因区域
                for line4 in open(Ref):#找出主要转录本选出注释区域
                    GeneRef=[]
                    if '%s\t'%gene[0] in line4:#获得主表达转录本
                        #print gene[0],line4.strip().split('\t')[13]
                        GeneRef.append(line4.strip().split('\t')[13])
                        #print gene,GeneRef
                    if len(GeneRef)>0:
                        ExonSite=[]
                        for line5 in open(gff):
                            if '%s"'%GeneRef[0] in line5:
                                if '\texon\t' in line5:
                                    #print gene,GeneRef,line5
                                    ExonSite.append(int(line5.split('\t')[3]))
                                    ExonSite.append(int(line5.split('\t')[4]))
                        ExonSite=list(set(ExonSite));ExonSite.sort()
                        #print ExonSite
                        #去除超出起始位点的PAS
                        while min(site)<min(ExonSite):#上游位点
                            site.remove(min(site))
                        while max(site)-10000>max(ExonSite):#下游超出10000bp的终止位点放弃
                            site.remove(max(site))
                        if len(site)>=2:#去除了PAS位点超出转录本起始位点后，数量还是2以上的进行分析
                            MinPAS=min(site);MaxPAS=max(site)
                            if min(ExonSite)<=MinPAS<=MaxPAS<=max(ExonSite):#全在范围内的只需要计算截取最近端PAS所在外显子及其后方的区域进行输出
                                #print 'All site in Range'#print 'ExonSite',ExonSite#print 'gene,MinPAS,MaxPAS',gene,MinPAS,MaxPAS
                                n=0
                                while n<len(ExonSite):
                                    #print ExonSite[n],ExonSite[n+1]#开始判定最近端PAS所在的exon进行输出
                                    if ExonSite[n+1]<MinPAS:#右端位点都小于minPAS的直接跳过
                                        pass
                                        n+=2
                                    elif ExonSite[n]<=MinPAS<=ExonSite[n+1]:
                                        #print '%s\tPacBio\tgene\t%s\t%s\t.\t+\t.\tgene "%s"\n'%(CHR[0],ExonSite[n],ExonSite[n+1],gene[0])
                                        myfile.write('%s\tPacBio\tgene\t%s\t%s\t.\t+\t.\tgene "%s"\n'%(CHR[0],ExonSite[n],ExonSite[n+1],gene[0]))
                                        n+=2
                                    elif MinPAS<=ExonSite[n]:
                                        #print '%s\tPacBio\tgene\t%s\t%s\t.\t+\t.\tgene "%s"\n'%(CHR[0],ExonSite[n],ExonSite[n+1],gene[0])
                                        myfile.write('%s\tPacBio\tgene\t%s\t%s\t.\t+\t.\tgene "%s"\n'%(CHR[0],ExonSite[n],ExonSite[n+1],gene[0]))
                                        n+=2
                                    else:
                                        print '?'*500
                                        n+=2
                            elif min(ExonSite)<=MinPAS and max(ExonSite)<MaxPAS:#最大位点超出范围，将最大值替换为MAXPAS进行判定
                                #print 'gene,MinPAS,MaxPAS',gene,MinPAS,MaxPAS
                                #print 'Before maxPAS',ExonSite
                                ExonSite.remove(max(ExonSite));ExonSite.append(MaxPAS)
                                ExonSite.sort()
                                #print 'After maxPAS',ExonSite
                                #print 'MaxPAS out of range'
                                n=0
                                while n<len(ExonSite):
                                    if ExonSite[n+1]<MinPAS:#右端位点都小于minPAS的直接跳过
                                        pass
                                        n+=2
                                    elif ExonSite[n]<=MinPAS<=ExonSite[n+1]:
                                        #print '%s\tPacBio\tgene\t%s\t%s\t.\t+\t.\tgene "%s"\n'%(CHR[0],ExonSite[n],ExonSite[n+1],gene[0])
                                        myfile.write('%s\tPacBio\tgene\t%s\t%s\t.\t+\t.\tgene "%s"\n'%(CHR[0],ExonSite[n],ExonSite[n+1],gene[0]))
                                        n+=2
                                    elif MinPAS<=ExonSite[n]:
                                        #print '%s\tPacBio\tgene\t%s\t%s\t.\t+\t.\tgene "%s"\n'%(CHR[0],ExonSite[n],ExonSite[n+1],gene[0])
                                        myfile.write('%s\tPacBio\tgene\t%s\t%s\t.\t+\t.\tgene "%s"\n'%(CHR[0],ExonSite[n],ExonSite[n+1],gene[0]))
                                        n+=2
                                    else:
                                        print '?'*500
                                        n+=2
                            elif MinPAS<min(ExonSite) and MaxPAS<=max(ExonSite):#左端位点超出范围，将最小值替换为MinPAS，无需比对直接输出
                                #print 'gene,MinPAS,MaxPAS',gene,MinPAS,MaxPAS
                                #print 'Before maxPAS',ExonSite
                                ExonSite.remove(min(ExonSite));ExonSite.append(MinPAS)
                                ExonSite.sort()
                                #print 'After maxPAS',ExonSite
                                n=0
                                while n<len(ExonSite):
                                    if ExonSite[n+1]<MinPAS:#右端位点都小于minPAS的直接跳过
                                        pass
                                        n+=2
                                    elif ExonSite[n]<=MinPAS<=ExonSite[n+1]:
                                        #print '%s\tPacBio\tgene\t%s\t%s\t.\t+\t.\tgene "%s"\n'%(CHR[0],ExonSite[n],ExonSite[n+1],gene[0])
                                        myfile.write('%s\tPacBio\tgene\t%s\t%s\t.\t+\t.\tgene "%s"\n'%(CHR[0],ExonSite[n],ExonSite[n+1],gene[0]))
                                        n+=2
                                    elif MinPAS<=ExonSite[n]:
                                        #print '%s\tPacBio\tgene\t%s\t%s\t.\t+\t.\tgene "%s"\n'%(CHR[0],ExonSite[n],ExonSite[n+1],gene[0])
                                        myfile.write('%s\tPacBio\tgene\t%s\t%s\t.\t+\t.\tgene "%s"\n'%(CHR[0],ExonSite[n],ExonSite[n+1],gene[0]))
                                        n+=2
                                    else:
                                        print '?'*500
                                        n+=2
                                #print 'MinPAS out of range'*100
                            elif MinPAS<min(ExonSite) and max(ExonSite)<MaxPAS:#左右两端都超出范围，最大最小位点都替换
                                ExonSite.remove(min(ExonSite));ExonSite.append(MinPAS); ExonSite.remove(max(ExonSite));ExonSite.append(MaxPAS)
                                ExonSite.sort()
                                #print 'After maxPAS',ExonSite
                                #print 'MaxPAS out of range'
                                n=0
                                while n<len(ExonSite):
                                    if ExonSite[n+1]<MinPAS:#右端位点都小于minPAS的直接跳过
                                        pass
                                        n+=2
                                    elif ExonSite[n]<=MinPAS<=ExonSite[n+1]:
                                        #print '%s\tPacBio\tgene\t%s\t%s\t.\t+\t.\tgene "%s"\n'%(CHR[0],ExonSite[n],ExonSite[n+1],gene[0])
                                        myfile.write('%s\tPacBio\tgene\t%s\t%s\t.\t+\t.\tgene "%s"\n'%(CHR[0],ExonSite[n],ExonSite[n+1],gene[0]))
                                        n+=2
                                    elif MinPAS<=ExonSite[n]:
                                        #print '%s\tPacBio\tgene\t%s\t%s\t.\t+\t.\tgene "%s"\n'%(CHR[0],ExonSite[n],ExonSite[n+1],gene[0])
                                        myfile.write('%s\tPacBio\tgene\t%s\t%s\t.\t+\t.\tgene "%s"\n'%(CHR[0],ExonSite[n],ExonSite[n+1],gene[0]))
                                        n+=2
                                    else:
                                        print '?'*500
                            else:
                                print 'What ?'*100
            if strand[0]=='-':
                #MinPAS=max(site);MaxPAS=min(site)
                #print gene,MinPAS,MaxPAS
                for JB in site:
                    myfile.write('%s\tPacBio\tapa\t%s\t%s\t.\t-\t.\tapa "%s.%s_%s"\n'%(CHR[0],JB,JB,site.index(JB),JB,gene[0]))
                #编写基因区域
                for line6 in open(Ref):
                    GeneRef=[]
                    if '%s\t'%gene[0] in line6:#获得主表达转录本
                        #print gene[0],line6.strip().split('\t')[13]
                        GeneRef.append(line6.strip().split('\t')[13])
                    if len(GeneRef)>0:
                        ExonSite=[]
                        for line7 in open(gff):
                            if '%s"'%GeneRef[0] in line7:
                                if '\texon\t' in line7:
                                        #print gene,GeneRef,line5
                                        ExonSite.append(int(line7.split('\t')[3]))
                                        ExonSite.append(int(line7.split('\t')[4]))
                        ExonSite=list(set(ExonSite));ExonSite.sort()
                        #print GeneRef,ExonSite
                        ExonSite.sort();ExonSite.reverse()
                        #去除超出起始位点的PAS
                        while max(site)>max(ExonSite):
                            site.remove(max(site))
                        while min(site)+10000<min(ExonSite):
                            site.remove(min(site))
                        if len(site)>=2:#去除了PAS位点超出转录本起始位点后，数量还是2以上的进行分析
                            MinPAS=max(site);MaxPAS=min(site)
                            #print GeneRef,ExonSite
                            if min(ExonSite)<=MaxPAS<=MinPAS<=max(ExonSite):#反向的统计
                                #print 'All site in Range'
                                n=0
                                while n<len(ExonSite):
                                    if ExonSite[n+1]>MinPAS:#左端位点都大于minPAS的直接跳过
                                        pass
                                        n+=2
                                    elif ExonSite[n+1]<=MinPAS<=ExonSite[n]:
                                        #print '%s\tPacBio\tgene\t%s\t%s\t.\t+\t.\tgene "%s"\n'%(CHR[0],ExonSite[n+1],ExonSite[n],gene[0])
                                        myfile.write('%s\tPacBio\tgene\t%s\t%s\t.\t-\t.\tgene "%s"\n'%(CHR[0],ExonSite[n+1],ExonSite[n],gene[0]))
                                        n+=2
                                    elif MinPAS>=ExonSite[n+1]:
                                        #print '%s\tPacBio\tgene\t%s\t%s\t.\t+\t.\tgene "%s"\n'%(CHR[0],ExonSite[n+1],ExonSite[n],gene[0])
                                        myfile.write('%s\tPacBio\tgene\t%s\t%s\t.\t-\t.\tgene "%s"\n'%(CHR[0],ExonSite[n+1],ExonSite[n],gene[0]))
                                        n+=2
                                    else:
                                        print '?'*500
                                        n+=2                  
                            elif max(ExonSite)>=MinPAS and min(ExonSite)>MaxPAS:#最大位点超出位点进行替换
                                ExonSite.remove(min(ExonSite));ExonSite.append(MaxPAS)
                                ExonSite.sort();ExonSite.reverse()
                                #print 'MaxPAS out of range'
                                n=0
                                while n<len(ExonSite):
                                    if ExonSite[n+1]>MinPAS:#左端位点都大于minPAS的直接跳过
                                        pass
                                        n+=2
                                    elif ExonSite[n+1]<=MinPAS<=ExonSite[n]:
                                        #print '%s\tPacBio\tgene\t%s\t%s\t.\t+\t.\tgene "%s"\n'%(CHR[0],ExonSite[n+1],ExonSite[n],gene[0])
                                        myfile.write('%s\tPacBio\tgene\t%s\t%s\t.\t-\t.\tgene "%s"\n'%(CHR[0],ExonSite[n+1],ExonSite[n],gene[0]))
                                        n+=2
                                    elif MinPAS>=ExonSite[n+1]:
                                        #print '%s\tPacBio\tgene\t%s\t%s\t.\t+\t.\tgene "%s"\n'%(CHR[0],ExonSite[n+1],ExonSite[n],gene[0])
                                        myfile.write('%s\tPacBio\tgene\t%s\t%s\t.\t-\t.\tgene "%s"\n'%(CHR[0],ExonSite[n+1],ExonSite[n],gene[0]))
                                        n+=2
                                    else:
                                        print '?'*500
                                        n+=2
                            elif MaxPAS>=min(ExonSite) and max(ExonSite)<MinPAS:
                                ExonSite.remove(max(ExonSite));ExonSite.append(MinPAS)
                                ExonSite.sort();ExonSite.reverse()
                                #print 'MinPAS out of range'*100
                                n=0
                                while n<len(ExonSite):
                                    if ExonSite[n+1]>MinPAS:#左端位点都大于minPAS的直接跳过
                                        pass
                                        n+=2
                                    elif ExonSite[n+1]<=MinPAS<=ExonSite[n]:
                                        #print '%s\tPacBio\tgene\t%s\t%s\t.\t+\t.\tgene "%s"\n'%(CHR[0],ExonSite[n+1],ExonSite[n],gene[0])
                                        myfile.write('%s\tPacBio\tgene\t%s\t%s\t.\t-\t.\tgene "%s"\n'%(CHR[0],ExonSite[n+1],ExonSite[n],gene[0]))
                                        n+=2
                                    elif MinPAS>=ExonSite[n+1]:
                                        #print '%s\tPacBio\tgene\t%s\t%s\t.\t+\t.\tgene "%s"\n'%(CHR[0],ExonSite[n+1],ExonSite[n],gene[0])
                                        myfile.write('%s\tPacBio\tgene\t%s\t%s\t.\t-\t.\tgene "%s"\n'%(CHR[0],ExonSite[n+1],ExonSite[n],gene[0]))
                                        n+=2
                                    else:
                                        print '?'*500
                                        n+=2
                            elif MaxPAS<min(ExonSite) and max(ExonSite)<MinPAS:
                                ExonSite.remove(min(ExonSite));ExonSite.append(MaxPAS);ExonSite.remove(max(ExonSite));ExonSite.append(MinPAS)
                                ExonSite.sort();ExonSite.reverse()
                                while n<len(ExonSite):
                                    if ExonSite[n+1]>MinPAS:#左端位点都大于minPAS的直接跳过
                                        pass
                                        n+=2
                                    elif ExonSite[n+1]<=MinPAS<=ExonSite[n]:
                                        #print '%s\tPacBio\tgene\t%s\t%s\t.\t+\t.\tgene "%s"\n'%(CHR[0],ExonSite[n+1],ExonSite[n],gene[0])
                                        myfile.write('%s\tPacBio\tgene\t%s\t%s\t.\t-\t.\tgene "%s"\n'%(CHR[0],ExonSite[n+1],ExonSite[n],gene[0]))
                                        n+=2
                                    elif MinPAS>=ExonSite[n+1]:
                                        #print '%s\tPacBio\tgene\t%s\t%s\t.\t+\t.\tgene "%s"\n'%(CHR[0],ExonSite[n+1],ExonSite[n],gene[0])
                                        myfile.write('%s\tPacBio\tgene\t%s\t%s\t.\t-\t.\tgene "%s"\n'%(CHR[0],ExonSite[n+1],ExonSite[n],gene[0]))
                                        n+=2
                                    else:
                                        print '?'*500
                                        n+=2
                                #print 'MinPAS and MAxPAS all out of range'*100
                                
                            else:
                                print 'What ?'*100

    myfile.close()

#roar2('Already-Fusion-PH1-fusion-version2.3.gtf.gff3.gtf-polyA-site-count.txt','All-PAS-Number-count.txt','Gene-T-for-Allgene-add-noCDS.txt','PH1-fusion-version4.2.gtf')
#去除基因区域外位点
def noOverAPA(x):
    myfile=open(r'roar3.gtf','w')
    myfile1=open(r'roar3-igv.gtf','w')
    geneid=[]
    genelist=[]#不知哪里出现了geneid重复，去除
    for line1 in open(x):
        if 'gene\t' in line1:
            if line1.split('"')[1] not in geneid:
                geneid.append(line1.split('"')[1])
    for i in geneid:
        print i
        Range=[]
        for line2 in open(x):
            if 'gene\t' in line2:
                if '%s"'%i in line2:
                    Range.append(int(line2.split('\t')[3]))
                    Range.append(int(line2.split('\t')[4]))
        for line3 in open(x):
            if 'apa\t' in line3:
                if '_%s'%i in line3:
                    if min(Range)<=int(line3.split('\t')[3])<=max(Range):#只要在范围内的
                        myfile.write(line3)
                        myfile1.write(line3)
            if 'gene\t' in line3:
                if '%s"'%i in line3:
                    if line3 not in genelist:
                        myfile.write(line3)
                        myfile1.write(line3.replace('gene\t','exon\t',1))
                        genelist.append(line3)
                            
    myfile.close()
    myfile1.close()
#noOverAPA('roar2.gtf')

#统计过滤后只剩一个位点的基因
def countAPA(x):
    myfile=open('roar-4.gtf','w')
    geneid=[]
    for line1 in open(x):
        if 'gene\t' in line1:
            if line1.split('"')[1] not in geneid:
                geneid.append(line1.split('"')[1])
    Single=[]
    for i in geneid:
         #print i
         APA=[]#统计APA数量
         for line3 in open(x):
            if 'apa\t' in line3:
                if '_%s'%i in line3:
                    APA.append(line3.split('"')[1])
         if len(APA)==1:
            Single.append(i)
    for line4 in open(x):
        if 'gene\t' in line4:
            if line4.split('"')[1] not in Single:
                myfile.write(line4)
        if 'apa\t' in line4:
            if line4.split('"')[1].split('_')[1] not in Single:
                myfile.write(line4)
    print len(Single)
    myfile.close()

#countAPA('roar3.gtf')

#检测有基因没为点的和有位点没基因的
def HavePASnoGene(x):
    Gene=[]
    for line1 in open(x):
        if 'gene\t' in line1:
            Gene.append(line1.split('gene "')[1].split('"')[0])
    Gene=list(set(Gene))
    #print len(Gene)
    for i in Gene:
        #print i
        APA=[]
        for line2 in open(x):
            if 'apa\t' in line2:
                if '_%s"'%i in line2:
                    APA.append(line2.split('apa "')[0])
        if len(APA)==0:
            print i,'Have no APA site'
HavePASnoGene('roar-4.gtf')
def HavePASnoGene2(x):
    APA=[]
    for line1 in open(x):
        if 'apa\t' in line1:
            APA.append(line1.split('_')[1].split('"')[0])
            #print line1.split('_')[1].split('"')[0]
    APA=list(set(APA))
    for i in APA:
        Gene=[]
        for line2 in open(x):
            if 'gene\t' in line2:
                if '%s"'%i in line2:
                    Gene.append(line2.split('gene "')[1].split('"')[0])
        if len(Gene)==0:
            print i,'Have no Gene site'

#HavePASnoGene2('roar-4.gtf')
