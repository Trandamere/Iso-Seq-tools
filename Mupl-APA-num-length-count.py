#coding:utf-8
#part7#将已有基因的转录本终止位点进行长度和数量统计
#逐个输出三代基因,
import numpy
import os
import sys

    
#生成拆分文件并分别运行程序
def splitFile(File,N,tracking,PASsite):
        Gene=[]
        for line1 in open(File):
                Gene.append( line1.split('gene_id "')[1].split('"')[0])
        Gene=list(set(Gene))
        Gene.sort()
        print len(Gene)#要拆分的总数
        print N#要拆分的份数
        Num=len(Gene)/N
        print len(Gene)/N#拆分的每一份
        #print len(Gene)-len(Gene)/N*N#不能正拆的余数
        list_of_groups = zip(*(iter(Gene),) *Num)#组合列表和差分数
        end_list = [list(i) for i in list_of_groups]
        count = len(Gene) %Num
        end_list.append(Gene[-count:]) if count !=0 else end_list
        Flag=1
        #根据endlist长度创建相应数量的文件，并根据文件编号进入所对应的拆分子列表，循环基因号打印文件
        #生成相应文件数量
        Fnum=len(end_list)
        print Fnum
        while Flag <=Fnum:
            myfile=open('F%s.gtf'%Flag,'w')
            #根据位置索引将向对应部分的拆分基因号遍历
            for Geneid in end_list[Flag-1]:
                for lineA in open(File):
                    if '%s"'%Geneid in lineA:
                        #print Geneid,lineA
                        myfile.write(lineA)
            Flag+=1
            myfile.close()
        
        #gtf拆分文件完成，接下来进行每个文件的分开现成运行
        Flag2=1
        Fnum2=len(end_list)
        while Flag2 <=Fnum2:
            print 'F%s.gtf'%Flag2
            os.system('nohup python flnc-chr.py F%s.gtf %s %s &'%(Flag2,tracking,PASsite))
            #print A
            Flag2+=1
            #os.system('python'%(x,x))
            #delNAB('FusionFG1G2-PH1-fusion-version2.3.gtf.gff3.gtf','PH1-flnc-mix12-hy12-sex12.tracking','FilterArich-flnc-mix12-hy12-sex12.gtf')

splitFile('Fusion-PH1-fusion-version2.3.gtf.gff3.gtf',14,'PH1-flnc-sex12.tracking','FilterArich-isoseq_flnc-All-Sex.gff3.gtf')
#splitFile(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
