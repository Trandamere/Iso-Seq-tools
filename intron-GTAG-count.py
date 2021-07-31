#!/usr/bin/env python3
# coding:utf-8
import argparse
import os
import sys
parser = argparse.ArgumentParser(description="功能：根据输入的gtf注释及基因组文件，生成内含子注释文件，提取内含子序列，提取内含子GTAG序列，进行GTAG统计。需要使用软件cufflinks")
parser.add_argument('-gtf', "--gtf", required=True,  help="gtf文件所在路径")
parser.add_argument('-genome', "--genome", required=True,  help="基因组序列所在路径")
parser.add_argument('-o', "--out", required=True,  help="数据输出文件夹路径，如果目录不存在则会自动创建")
args = parser.parse_args()
out_dir = os.path.dirname(args.out)
if out_dir and not os.path.exists(out_dir):
    os.makedirs(out_dir)

#首先进行内含子文件生成
def intronGff(gff,genome,out):
    print('Creating intron transcript gtf')
    myfile=open('%stran-intron.gff'%out,'w')
    Transcript={}
    for line in open(gff):
        if 'exon\t' in line:
            chro=line.split('\t')[0]
            lian=line.split('\t')[6]
            gene=line.split('"')[1]
            T=line.split('"')[3]
            #print Transcript
            Site1=int(line.split('\t')[3])
            Site2=int(line.split('\t')[4])           
            #此处进行判定，键存在的话，替换值的第4个，进行扩充
            if T in Transcript.keys():##进行位点累积替换算法
                Site=sorted([i for i in Transcript[T][3]])
                Site.append(Site1);Site.append(Site2)
                Transcript[T]=[chro,lian,gene,Site]
            else:
                Transcript[T]=[chro,lian,gene,[Site1,Site2]]
    for key,values in Transcript.items():
        #print key,values
        chro=values[0];Exon=values[3];del Exon[0];del Exon[-1];pop=key;gene=values[2];lian=values[1]
        n=0
        while n+1<=len(Exon):
            intron_start=Exon[n]+1;intron_stop=Exon[n+1]-1
            length =intron_stop-intron_start+1
            n+=2
            myfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tgene_id "%s"; transcript_id "%s-%s-%s-len.%s";\n' %(chro,'luping','exon',intron_start,intron_stop,'.',lian,'.',gene,pop,intron_start,intron_stop,length))
    myfile.close()
    print('Done')
    ##############################################################################33
    myfile1=open('%sgene-intron.gff'%out,'w')
    print('Creating gene intron transcript gtf')
    Gene=[]#修改为字典，将内含子添加入基因为键，id为值的字典，如果有重复过滤，不重复则添加，完成后。遍历基因的键值，输出基因的内含子
    for line1 in open('%stran-intron.gff'%out):
        part1=line1.split('transcript_id "')[0];ID=line1.split('"')[1];part2=line1.split('transcript_id "')[1].replace(line1.split('transcript_id "')[1].split('-')[0],ID)
        #print line1
        GeneIntron= '%stranscript_id "%s'%(part1,part2)
        Gene.append(GeneIntron)
    #print len(Gene)
    Gene=sorted(list(set(Gene)))
    #print len(Gene)
    for i in Gene:
        myfile1.write(i)
    myfile1.close()
    print('Done')
intronGff(args.gtf,args.genome,args.out)
