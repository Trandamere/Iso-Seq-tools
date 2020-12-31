#coding:utf-8
#part1
import sys
import numpy
def delNAB(gff,tracking,PASsite):
        myfile1=open('%s-PAS-number-count.txt'%gff,'w')
        myfile2=open('%s-PAS-Length-count.txt'%gff,'w')
        myfile1.write('ID\tPAS\tNumber\n')
        myfile2.write('ID\tPASlist\tLength\n')
        Same=['=','c','e','j','k','o']
        #g根据基因号码，逐个输出基因号，调取tracking比对结果，调取覆盖转录本，提取终止位点为PAS位点，进行数量统计与距离计算
        Gene=[]
        for line1 in open(gff):
                Gene.append( line1.split('gene_id "')[1].split('"')[0])
        Gene=list(set(Gene))
        Gene.sort()
        for i in Gene:
                print i
                Fastanumber=[]
                for line2 in open(tracking):
                        if '%s|'%i in line2:
                                if line2.split('\t')[3] in Same:
                                        Fastanumber.append(line2.split('q1:')[1].split('|')[1])
                                        
                if len(Fastanumber)>1:#有一个以上的初始序列比对，才能进行统计
                        PAScluster=[]
                        for j in Fastanumber:#遍历序列名称，进行位点提取
                                
                                for line3 in open(PASsite):
                                        if '%s"'%j in line3:
                                                #print i,j,line3
                                                if 'transcript\t' in line3:
                                                        if line3.split('\t')[6]=='+':
                                                                PAScluster.append(int(line3.split('\t')[4]))
                                                        if line3.split('\t')[6]=='-':
                                                                PAScluster.append(int(line3.split('\t')[3]))
                        PAScluster2=list(set(PAScluster))#第一部分#输出位点统计及支持的数量文件
                        PAScluster2.sort()
                        for k in PAScluster2:
                                PAScount=PAScluster.count(k)
                                #print i,'%s-%s'%(i,k),'Num',PAScount#文件数量统计
                                Genesite='%s-%s'%(i,k)
                                myfile1.write('%s\t%s\t%s\n'%(i,Genesite,PAScount))
                                #print i,j,Fastanumber,PAScluster,k,PAScount,'Site Count'*10
                        #print i,j,Fastanumber,PAScluster
                        PAScluster3=list(set(PAScluster))#进行错位相减，统计长度
                        PAScluster4=list(set(PAScluster))
                        if len(PAScluster3)>1:#位点去重复后保证有两个位点#第二部分#输出临近位点
                                PAScluster3.sort();PAScluster4.sort()
                                del PAScluster3[-1];del PAScluster4[0]
                                #print PAScluster3
                                #print PAScluster4
                                PAScluster3=numpy.array(PAScluster3)
                                PAScluster4=numpy.array(PAScluster4)
                                L=PAScluster4-PAScluster3
                                for I in L:
                                        #print i,PAScluster2,'Length',I
                                        myfile2.write('%s\t%s\t%s\n'%(i,PAScluster2,I))
        myfile1.close()
        myfile2.close()
                        
#delNAB('Fusion-PH1-fusion-version2.3.gtf.gff3.gtf','PH1-ISO.tracking','FilterArich-hlq-flag-fix.gff3.gtf')
#delNAB('FusionFG2G-PH1-fusion-version2.3.gtf.gff3.gtf','PH1-ISO-mix12-hy12-sex12.tracking','FilterArich-mix12-hy12-sex12.gtf')
#delNAB('FusionFG2G-PH1-fusion-version2.3.gtf.gff3.gtf','PH1-flnc-mix12-hy12-sex12.tracking','FilterArich-flnc-mix12-hy12-sex12.gtf')
delNAB(sys.argv[1],sys.argv[2],sys.argv[3])
