#coding:utf-8
#2020.8.11
#luping
#查找APA
import os
import numpy
#根据三代数据末端为polyA的特征，进行APA分析，
#原理为将每个基因号下属转录本终止位点收集，从近到远，以最近的位点为起点，以6bp
#为范围划分为同一个polyA位点
def APA(Countfile,GFF,GFForg):
    myfile1=open(r'%s-polyA-site-count.txt'%GFF,'w')
    myfile1.write('Geneid\tCHR\tStrand\tPolyA-number\tLocation\tSite-cluster\n')
    Gene=[];
    for line1 in open(Countfile):
        if 'FG' in line1:
            Gene.append( line1.split('\t')[0])
    Gene=list(set(Gene))
    Gene.sort()
    #print len(Gene),Gene[0]
    for FGid in Gene:
        print FGid
        tranid=[];strand=[];polyAsite=[];CHR=[]
        for line2 in open(GFF):#收集链方向和染色体
            if '%s"'%FGid in line2:
                    if 'transcript\t' in line2:
                        strand.append(line2.split('\t')[6])
                        CHR.append(line2.split('\t')[0])
        #print strand[0],CHR[0]
        for line3 in open(Countfile):#收集基因位点信息：
            if '%s\t'%FGid in line3:
                Sites=line3.split('\t')[1].split('[')[1].split(']')[0].split(', ')
                for i in Sites:
                    polyAsite.append(int(i))
        polyAsite=list(set(polyAsite))
        polyAsite.sort()
        for j in polyAsite:
            #print j,strand[0],CHR[0]
            T=[]
            for line4 in open(GFForg):
                if 'transcript\t' in line4:
                    if '\t%s\t'%j in line4:
                        if line4.split('\t')[6]==strand[0]:
                            if line4.split('\t')[0]==CHR[0]:
                                T.append( line4.split('transcript_id "')[1].split('"')[0])
            T2="|".join(T)
            #print T
            #print T2
            tranid.append(T2)
        #print FGid,tranid,polyAsite
        tranApoly=[]#组合转录本和位点
        while len(tranid)>0:#收集基因下属转录本及其polyA位点
            tranidpop=tranid.pop();polyAsitepop=polyAsite.pop()
            T=[tranidpop,polyAsitepop]
            tranApoly.append(T)
        #print tranApoly
        if len(strand)>0:
            if strand[0]=='+':#根据链方向调整最近的PAS位点排序#反向，后续倒装
                tranApoly2=sorted(tranApoly,key=lambda x:(-x[1]))
            else:
                tranApoly2=sorted(tranApoly,key=lambda x:(x[1]))
            #tranApoly2=sorted(tranApoly,key=lambda x:(x[1]))
            #print tranApoly2
                
            if len(tranApoly2)==0:
                pass
            elif len(tranApoly2)==1:
                #print '%s\t%s\t%s\t1\t%s\t%s\n'%(FGid,CHR[0],strand[0],tranApoly2[0][1],tranApoly2)
                myfile1.write('%s\t%s\t%s\t1\t%s\t%s\n'%(FGid,CHR[0],strand[0],tranApoly2[0][1],tranApoly2))
            else:
                ID=[];number=[];SITE=[];Location=[]#将各个位点信息集合一起，方便统计总的位点信息
                n=1
                while len(tranApoly2)>0:
                    if strand[0]=='+':
                        number.append(n)#将位点数添加
                        pop=tranApoly2.pop()#原列表删除倒数第一个
                        same=[]
                        same.append(pop)#取出倒数第一个
                        Range=[]
                        Range.append(pop[1])
                        #site1MAX=same[-1][1]+20#最大值为去除的第一个
                        site1MAX=same[-1][1]+9#最大值为去除的第一个
                        #print 'same',same,'pop',pop,'tranApoly2',tranApoly2,'site1MAX',site1MAX
                        while len(tranApoly2)>0:#当剩余转录本还有
                            if site1MAX>=tranApoly2[-1][1]:#如果剩余转录本下一个临近位点的数值在20bp范围内
                                pop2=tranApoly2.pop()#取出这个位点
                                same.append(pop2)#添加到同一位点中
                                Range.append(pop2[1])#添加变化范围
                                site1MAX=site1MAX-same[-2][1]+same[-1][1]#换成下一个近端位点
                                #print 'same',same,'pop',pop,'tranApoly2',tranApoly2,'site1MAX',site1MAX,'!'*50
                            else:
                                break#没有位点则跳出循环
                        
                        SITE.append(same)
                        Location.append(Range)
                        n+=1

                    if strand[0]=='-':
                        number.append(n)#将位点数添加
                        pop=tranApoly2.pop()#原列表删除倒数第一个
                        same=[]
                        same.append(pop)#取出倒数第一个
                        Range=[]
                        Range.append(pop[1])
                        #site1MAX=same[-1][1]-20#最大值为去除的第一个#因为是反向，向前延伸20bp
                        site1MAX=same[-1][1]-9#最大值为去除的第一个#因为是反向，向前延伸6bp
                        #print 'same',same,'pop',pop,'tranApoly2',tranApoly2,'site1MAX',site1MAX
                        while len(tranApoly2)>0:#当剩余转录本还有
                            if site1MAX<=tranApoly2[-1][1]:#如果剩余转录本下一个临近位点的数值在20bp范围内#因为反向改为小于
                                pop2=tranApoly2.pop()#取出这个位点
                                same.append(pop2)#添加到同一位点中
                                Range.append(pop2[1])#添加变化范围
                                site1MAX=site1MAX-same[-2][1]+same[-1][1]#换成下一个近端位点
                                #print 'same',same,'pop',pop,'tranApoly2',tranApoly2,'site1MAX',site1MAX,'!'*50
                            else:
                                break#没有位点则跳出循环
                        
                        SITE.append(same)
                        Location.append(Range)
                        n+=1

                #print FGid,max(number),Location,SITE
                #print '%s\t%s\t%s\t%s\t%s\t%s\n'%(FGid,CHR[0],strand[0],max(number),Location,SITE)
                myfile1.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(FGid,CHR[0],strand[0],max(number),Location,SITE))
                                     
    myfile1.close()

#ALL时期运行
#cluster距离还是6不更改，更换A富集去除序列
#APA('All-PAS-Length-count.txt','Fusion-PH1-fusion-version2.3.gtf.gff3.gtf','FilterArich-mix12-hy12-sex12.gtf')
    

#cluster距离改成9#中途中断，加载计算完的
#APA('All-PAS-Length-count.txt','Fusion-PH1-fusion-version2.3.gtf.gff3.gtf','FilterArich-flnc-mix12-hy12-sex12.gtf')
def APA2(Already,Countfile,GFF,GFForg):
    myfile1=open(r'%s-polyA-site-count.txt'%GFF,'w')
    myfile1.write('Geneid\tCHR\tStrand\tPolyA-number\tLocation\tSite-cluster\n')
    AlreadyGene=[]
    for lineAlready in open(Already):
        AlreadyGene.append( lineAlready.split('\t')[0])
    Gene=[];
    for line1 in open(Countfile):
        if 'FG' in line1:
            Gene.append( line1.split('\t')[0])
    Gene=list(set(Gene))
    Gene.sort()
    #print len(Gene),Gene[0]
    for FGid in Gene:
        if FGid not in AlreadyGene:
            print FGid
            tranid=[];strand=[];polyAsite=[];CHR=[]
            for line2 in open(GFF):#收集链方向和染色体
                if '%s"'%FGid in line2:
                        if 'transcript\t' in line2:
                            strand.append(line2.split('\t')[6])
                            CHR.append(line2.split('\t')[0])
            #print strand[0],CHR[0]
            for line3 in open(Countfile):#收集基因位点信息：
                if '%s\t'%FGid in line3:
                    Sites=line3.split('\t')[1].split('[')[1].split(']')[0].split(', ')
                    for i in Sites:
                        polyAsite.append(int(i))
            polyAsite=list(set(polyAsite))
            polyAsite.sort()
            for j in polyAsite:
                #print j,strand[0],CHR[0]
                T=[]
                for line4 in open(GFForg):
                    if 'transcript\t' in line4:
                        if '\t%s\t'%j in line4:
                            if line4.split('\t')[6]==strand[0]:
                                if line4.split('\t')[0]==CHR[0]:
                                    T.append( line4.split('transcript_id "')[1].split('"')[0])
                T2="|".join(T)
                #print T
                #print T2
                tranid.append(T2)
            #print FGid,tranid,polyAsite
            tranApoly=[]#组合转录本和位点
            while len(tranid)>0:#收集基因下属转录本及其polyA位点
                tranidpop=tranid.pop();polyAsitepop=polyAsite.pop()
                T=[tranidpop,polyAsitepop]
                tranApoly.append(T)
            #print tranApoly
            if len(strand)>0:
                if strand[0]=='+':#根据链方向调整最近的PAS位点排序#反向，后续倒装
                    tranApoly2=sorted(tranApoly,key=lambda x:(-x[1]))
                else:
                    tranApoly2=sorted(tranApoly,key=lambda x:(x[1]))
                #tranApoly2=sorted(tranApoly,key=lambda x:(x[1]))
                #print tranApoly2
                    
                if len(tranApoly2)==0:
                    pass
                elif len(tranApoly2)==1:
                    #print '%s\t%s\t%s\t1\t%s\t%s\n'%(FGid,CHR[0],strand[0],tranApoly2[0][1],tranApoly2)
                    myfile1.write('%s\t%s\t%s\t1\t%s\t%s\n'%(FGid,CHR[0],strand[0],tranApoly2[0][1],tranApoly2))
                else:
                    ID=[];number=[];SITE=[];Location=[]#将各个位点信息集合一起，方便统计总的位点信息
                    n=1
                    while len(tranApoly2)>0:
                        if strand[0]=='+':
                            number.append(n)#将位点数添加
                            pop=tranApoly2.pop()#原列表删除倒数第一个
                            same=[]
                            same.append(pop)#取出倒数第一个
                            Range=[]
                            Range.append(pop[1])
                            #site1MAX=same[-1][1]+20#最大值为去除的第一个
                            site1MAX=same[-1][1]+9#最大值为去除的第一个
                            #print 'same',same,'pop',pop,'tranApoly2',tranApoly2,'site1MAX',site1MAX
                            while len(tranApoly2)>0:#当剩余转录本还有
                                if site1MAX>=tranApoly2[-1][1]:#如果剩余转录本下一个临近位点的数值在20bp范围内
                                    pop2=tranApoly2.pop()#取出这个位点
                                    same.append(pop2)#添加到同一位点中
                                    Range.append(pop2[1])#添加变化范围
                                    site1MAX=site1MAX-same[-2][1]+same[-1][1]#换成下一个近端位点
                                    #print 'same',same,'pop',pop,'tranApoly2',tranApoly2,'site1MAX',site1MAX,'!'*50
                                else:
                                    break#没有位点则跳出循环
                            
                            SITE.append(same)
                            Location.append(Range)
                            n+=1

                        if strand[0]=='-':
                            number.append(n)#将位点数添加
                            pop=tranApoly2.pop()#原列表删除倒数第一个
                            same=[]
                            same.append(pop)#取出倒数第一个
                            Range=[]
                            Range.append(pop[1])
                            #site1MAX=same[-1][1]-20#最大值为去除的第一个#因为是反向，向前延伸20bp
                            site1MAX=same[-1][1]-9#最大值为去除的第一个#因为是反向，向前延伸6bp
                            #print 'same',same,'pop',pop,'tranApoly2',tranApoly2,'site1MAX',site1MAX
                            while len(tranApoly2)>0:#当剩余转录本还有
                                if site1MAX<=tranApoly2[-1][1]:#如果剩余转录本下一个临近位点的数值在20bp范围内#因为反向改为小于
                                    pop2=tranApoly2.pop()#取出这个位点
                                    same.append(pop2)#添加到同一位点中
                                    Range.append(pop2[1])#添加变化范围
                                    site1MAX=site1MAX-same[-2][1]+same[-1][1]#换成下一个近端位点
                                    #print 'same',same,'pop',pop,'tranApoly2',tranApoly2,'site1MAX',site1MAX,'!'*50
                                else:
                                    break#没有位点则跳出循环
                            
                            SITE.append(same)
                            Location.append(Range)
                            n+=1

                    #print FGid,max(number),Location,SITE
                    #print '%s\t%s\t%s\t%s\t%s\t%s\n'%(FGid,CHR[0],strand[0],max(number),Location,SITE)
                    myfile1.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(FGid,CHR[0],strand[0],max(number),Location,SITE))
                                     
    myfile1.close()
#APA2('Already-Fusion-PH1-fusion-version2.3.gtf.gff3.gtf-polyA-site-count.txt','All-PAS-Length-count.txt','Fusion-PH1-fusion-version2.3.gtf.gff3.gtf','FilterArich-flnc-mix12-hy12-sex12.gtf')
APA2('Already-Fusion-PH1-fusion-version2.3.gtf.gff3.gtf-polyA-site-count.txt','delMt-All-PAS-Length-count.txt','delMt-Fusion-PH1-fusion-version2.3.gtf.gff3.gtf','delMt-FilterArich-flnc-mix12-hy12-sex12.gtf')
