#coding:utf-8
#2020.7.5
#将logMm值分成小于等-1，-1和1之间，大于等1
import sys
import numpy
def delNAB(tracking,PASsite):
        myfile1=open('%s-PAS-number-count.txt'%tracking,'w')
        myfile2=open('%s-PAS-Length-count.txt'%tracking,'w')
        myfile1.write('ID\tPAS\tNumber\n')
        myfile2.write('ID\tPASlist\tLength\n')
        Same=['=','c','e','j','k','o','m','n']
        #g根据基因号码，逐个输出基因号，调取tracking比对结果，调取覆盖转录本，提取终止位点为PAS位点，进行数量统计与距离计算
        #建立基因与序列关系的字典
        GeneFasta={}
        print('开始收集基因与css关系')
        for line2 in open(tracking):
                if line2.split('\t')[3] in Same:
                        PBID = line2.split('|')[0].split('\t')[2];CCSID = line2.split('q1:')[1].split('|')[1]
                        #将ID和css名字加入字典
                        if PBID not in GeneFasta.keys():
                                GeneFasta[PBID] = [CCSID]
                        else:
                                #如果有PBID，将css名字，加入字典中
                                Va = GeneFasta[PBID]
                                Va.append(CCSID)
                                GeneFasta[PBID] = Va
        # print(len(GeneFasta))
        print('收集完毕')
        # #建立序列名称和位点关系的字典
        FastaGff = {}#
        print('开始收集css与PAS位点关系')
        for line3 in open(PASsite):
                if 'transcript\t' in line3:
                        ID = line3.split('"')[3]
                        if line3.split('\t')[6] == '+':#每条转录本只有一个位点，不涉及多位点替换
                                FastaGff[ID] = int(line3.split('\t')[4])
                        if line3.split('\t')[6] == '-':  # 每条转录本只有一个位点，不涉及多位点替换
                                FastaGff[ID] = int(line3.split('\t')[3])
        # #所有CSS转录本位点末端信息字典记录
        # print('CSS转录本长度',len(FastaGff))
        print('收集完毕')
        # #现在有了PB号对应css名字的字典：#css名字对应位点的字典#现在需要建立一个新的字典，使用css作为中间链接部件，建立PB与位点的字典
        print('开始关联PB与PAS位点关系')
        GeneGff = {}#存放基因，PAS位点关系
        for key1,value1 in GeneFasta.items():
                PASlist = []
                for value2 in value1:#value2为css名称
                        # print(key1,value2,FastaGff[value2])
                        PASlist.append(FastaGff[value2])
                GeneGff[key1] = PASlist
        print('收集完毕')
        print('开始计算length和count文件')
        for key3,value3 in GeneGff.items():
                # print(key3,value3)
                value4 = list(set(value3))
                # print(key3, value4)
                for CoutN in value4:
                        PAScount = value3.count(CoutN)
                        #print(key3,value3,value4,CoutN,PAScount)
                        Genesite = '%s-%s' % (key3, CoutN)
                        myfile1.write('%s\t%s\t%s\n' % (key3, Genesite, PAScount))
                #将文件排序
                PAScluster2 = value4
                PAScluster2.sort()
                PAScluster3 = list(set(PAScluster2))  # 进行错位相减，统计长度
                PAScluster4 = list(set(PAScluster2))
                if len(PAScluster2)>1:
                        PAScluster3.sort();
                        PAScluster4.sort()
                        del PAScluster3[-1];
                        del PAScluster4[0]
                        PAScluster3 = numpy.array(PAScluster3)
                        PAScluster4 = numpy.array(PAScluster4)
                        L = PAScluster4 - PAScluster3
                        for I in L:
                                # print(key3,PAScluster2,I)
                                myfile2.write('%s\t%s\t%s\n' % (key3,PAScluster2,I))
        print('计算完成')
        myfile1.close()
        myfile2.close()
#APAeffectCDS('Fusion-PH1-fusion-version2.3.gtf.gff3.gtf','PH1-flnc-mix12-hy12-sex12.tracking','FilterArich-flnc-mix12-hy12-sex12.gtf')
delNAB(sys.argv[1],sys.argv[2])

