#coding:utf-8
import numpy
import time
import argparse
import os
import sys
parser = argparse.ArgumentParser(description="功能：根据emboss预测cds结果，gtf注释文件和intron注释文件输出匹配的cds注释文件。需要使用软件cufflinks")
parser.add_argument('-gtf', "--gtf", required=True,  help="gtf文件所在路径")
parser.add_argument('-cds', "--cds", required=True,  help="emboss 预测cds.fasta路径")
parser.add_argument('-intron', "--intron", required=True,  help="转录本内含子注释")
parser.add_argument('-o', "--out", required=True,  help="数据输出文件夹路径，如果目录不存在则会自动创建")
args = parser.parse_args()
out_dir = os.path.dirname(args.out)
if out_dir and not os.path.exists(out_dir):
    os.makedirs(out_dir)

def FGid(cds,gff,Intron,out):#根据cds提供的位置，计算在转录本的起始位置，
    myfile=open('%smax-CDS.gtf'%out,'w')
    Tsite={}#要exon位置信息和链方向和染色体编号进行计算
    CDSsite={}#只要起点和长度
    Intronsite={}
    #一次读取随后直接在里面全部计算，先读取cds和gtf信息

    for line1 in open(cds):
        if '>' in line1:
            CDSid=line1.split('_')[0].split('>')[1]
            CDSstart=int(line1.split(' - ')[0].split('[')[1])
            CDSstop=int(line1.split(' - ')[1].split(']')[0])
            CDSL=abs(int(line1.split(' - ')[0].split('[')[1])-int(line1.split(' - ')[1].split(']')[0]))
            CDSsite[CDSid]=[CDSstart,CDSstop,CDSL]
            #print CDSid,CDSstart,CDSL
            #print line1
    for line2 in open(gff):
        if 'exon\t' in line2:
            Tid=line2.split('transcript_id "')[1].split('"')[0]
            exon1=(int(line2.split('\t')[3]))
            exon2=(int(line2.split('\t')[4]))
            strand=(line2.split('\t')[6])
            CHR=(line2.split('\t')[0])
            Gene=(line2.split('"')[1])
            if Tid in Tsite.keys():
                Exon2=[i for i in Tsite[Tid][3]]
                Exon2.append(exon1);Exon2.append(exon2)
                Exon2.sort()
                Tsite[Tid]=[strand,CHR,Gene,Exon2]
            else:
                Exon=[exon1,exon2]
                Exon.sort()
                Tsite[Tid]=[strand,CHR,Gene,Exon]

    for line3 in open(Intron):
        Intronid=line3.split('transcript_id "')[1].split('-')[0]
        IntronSite=[int(line3.split('\t')[3]),int(line3.split('\t')[4])]
        if Intronid in Intronsite.keys():
            IntronSite2=[i for i in Intronsite[Intronid]]
            IntronSite2.append(IntronSite)
            Intronsite[Intronid]=IntronSite2
        else:
            Intronsite[Intronid]=[IntronSite]

    for key,value in Tsite.items():
        #print key,value
        if len(value[3])>2:#说明有内含子，进行计算
            if key in CDSsite.keys():#如果有结果
                #pass#随后进行计算
                #可以计算前后的内含子，计算起始或终止位点是否和剪切位点在同一个上
                #起始或终止位点在exon的起始或终止上需要添加这个点到CDS的位点上运行
                #可以一次性计算
                #print key,value,CDSsite[key],'intron number%s'%(len(value[3])/2-1)
                #查看起始位置是否与exon边界重合
                Already=[]
                if value[0]=='+':
                    TcdsStart=min(value[3])+CDSsite[key][0]-1
                    TcdsStop=min(value[3])+CDSsite[key][1]-1
                    #print(key,value[3])
                    #print('TcdsStart,TcdsStop1',TcdsStart,TcdsStop)
                    #print(key,value[3],CDSsite[key],Intronsite[key])
                    #接下来计算起点和终点位置累积的内含子长度
                    #############################################此处涉及到重新加长度的内含子怎么及时的和新的值比对
                    #此处做个while循环，内含子数量和内含子累计长度当场计算
                    n=0
                    while n<len(Intronsite[key]):
                        for I in Intronsite[key]:
                            if TcdsStart<I[0]:#cds起始位点在前内含子在后面，不用考虑内含子累加
                                pass
                                n+=1
                            elif TcdsStart>I[0] or I[1]:#cds起始位点大于内含子起始位点，把内含子长度累积加上去
                                Slength=abs(I[1]-I[0])+1
                                TcdsStart+=Slength
                                n+=1
                            elif TcdsStart==I[0]:#起始位点在exon点上，不算小于的内含子
                                pass
                                n+=1
                            else:
                                print(key,'what1?'*100)#除了大就是小，不应该有这部分
                                n+=1
                    m=0#同理起始位点，同样计算终止位点的距离
                    while m<len(Intronsite[key]):
                        for I in Intronsite[key]:
                            if TcdsStop<I[0]:#cds起始位点在前内含子在后面，不用考虑内含子累加
                                pass
                                m+=1
                            elif TcdsStop>I[0] or I[1]:#cds起始位点大于内含子起始位点，把内含子长度累积加上去
                                Slength=abs(I[1]-I[0])+1
                                TcdsStop+=Slength
                                m+=1
                            elif TcdsStop==I[0]:#终止位点在exon的内含子上，不算小于的内含子
                                pass
                                m+=1
                            else:
                                print(key,'what2?'*100)#除了大就是小，不应该有这部分
                                m+=1
                    #TcdsStop+=3#emboss未发表终止位点修正
                    #print('TcdsStart,TcdsStop2',TcdsStart,TcdsStop)
                    CDSexonSite=[]#存放在CDS范围内的exon
                    ##接下来把exon中归属CDS的全部放倒一个列表住，排列然后输出gtf
                    for II in value[3]:
                        if TcdsStart<II<TcdsStop:
                            CDSexonSite.append(II)
                    #此处起始位点在偶数位置时进行添加
                    if TcdsStart in (value[3]):#如果起始位点在exon的某个位点，则进行gtf需求式补充
                        if (value[3].index(TcdsStart))% 2 != 0:#添加索引修正后，是奇数则在偶数位置
                            CDSexonSite.append(TcdsStart)
                    #终止位点在技术位置时才进行添加
                    if TcdsStop in (value[3]):#如果终止位点在exon的某个位点，贼进行gtf需求时补充
                        if (value[3].index(TcdsStop))% 2 == 0:#添加索引修正后，是奇数则在偶数位置
                            CDSexonSite.append(TcdsStop+3)
                    CDSexonSite.append(TcdsStart);CDSexonSite.append(TcdsStop+3)#添加起点和终点
                    CDSexonSite.sort()
                    if len(CDSexonSite)% 2 == 0:
                        #pass#print len(CDSexonSite),CDSexonSite
                        N=0
                        while N<len(CDSexonSite):
                            #print CDSexonSite[N],CDSexonSite[N+1]
                            #print('%s\tPacBio\tCDS\t%s\t%s\t.\t%s\t.\tgene_id "%s"; transcript_id "%s";\n'%(value[1], CDSexonSite[N],CDSexonSite[N+1],value[0],value[2],key))
                            myfile.write('%s\tPacBio\tCDS\t%s\t%s\t.\t%s\t.\tgene_id "%s"; transcript_id "%s";\n'%(value[1], CDSexonSite[N],CDSexonSite[N+1],value[0],value[2],key))
                            N+=2
                    else:#位点不是成对的
                        print len(CDSexonSite),CDSexonSite
                        print(key,'?1'*1000)
                if value[0]=='-':
                    #print(key)#,CDSsite[key])
                    TcdsStart=max(value[3])-CDSsite[key][0]+1#负链的cds起始为最大exon减去cds起始的长度
                    TcdsStop=max(value[3])-CDSsite[key][1]+1#负链的cds终止为最大exon减去cds终止的长度
                    #print(TcdsStart,TcdsStop-3,value[3])
                    w=0
                    Intronsite[key].reverse()###此处需要累计计算内含子长度再进行比对，所以内含子输出顺序是关键，负链的内含子输出顺序是相反的
                    while w<len(Intronsite[key]):
                        for I in Intronsite[key]:
                            if I[1]<TcdsStart:#cds起始位点在前内含子在后面，不用考虑内含子累加
                                pass
                                w+=1
                            elif TcdsStart<I[0] or I[1]:#cds起始位点大于内含子起始位点，把内含子长度累积加上去
                                Slength=abs(I[1]-I[0])+1
                                TcdsStart-=Slength
                                w+=1
                            elif TcdsStart==I[1]:#起始位点在exon点上，不算小于的内含子
                                pass
                                w+=1
                            else:
                                print(key,'what3?'*100)#除了大就是小，不应该有这部分
                                w+=1
                    #print(TcdsStart)

                    v=0#同理起始位点，同样计算终止位点的距离
                    while v<len(Intronsite[key]):
                        for I in Intronsite[key]:
                            if I[1]<TcdsStop:#cds起始位点在前内含子在后面，不用考虑内含子累加
                                pass
                                v+=1
                            elif TcdsStop<I[1] or I[0]:#cds起始位点大于内含子起始位点，把内含子长度累积加上去
                                Slength=abs(I[1]-I[0])+1
                                TcdsStop-=Slength
                                v+=1
                            elif TcdsStop==I[1]:#起始位点在exon点上，不算小于的内含子
                                pass
                                v+=1
                            else:
                                print(key,'what4?'*100)#除了大就是小，不应该有这部分
                                v+=1
                    #print(TcdsStart,TcdsStop-3,value[3])
                    CDSexonSite=[]#存放在CDS范围内的exon
                    ##接下来把exon中归属CDS的全部放倒一个列表住，排列然后输出gtf
                    for III in value[3]:
                        if TcdsStop<III<TcdsStart:
                            CDSexonSite.append(III)
                    #print('CDSexonSite',CDSexonSite)
                    #此处起始位点在偶数位置时进行添加
                    #此处起始位点在偶数位置时进行添加
                    if TcdsStart in (value[3]):#如果起始位点在exon的某个位点，则进行gtf需求式补充#因为在负链，需要修正正向索引带来的影响
                        if (value[3].index(TcdsStart))% 2 == 0:#添加索引修正后，是偶数则在偶数位置
                            CDSexonSite.append(TcdsStart)
                    #终止位点在技术位置时才进行添加
                    if TcdsStop in (value[3]):#如果终止位点在exon的某个位点，则进行gtf需求时补充
                        if (value[3].index(TcdsStop))% 2 != 0:#添加索引修正后，是奇数则在偶数位置
                            CDSexonSite.append(TcdsStop-3)
                    CDSexonSite.append(TcdsStart);CDSexonSite.append(TcdsStop-3)#添加起点和终点
                    CDSexonSite.sort()
                    #print('CDSexonSite',CDSexonSite)
                    if len(CDSexonSite)% 2 == 0:
                        #pass#print len(CDSexonSite),CDSexonSite
                        N=0
                        while N<len(CDSexonSite):
                            #print CDSexonSite[N],CDSexonSite[N+1]
                            #print('%s\tPacBio\tCDS\t%s\t%s\t.\t%s\t.\tgene_id "%s"; transcript_id "%s";\n'%(value[1], CDSexonSite[N],CDSexonSite[N+1],value[0],value[2],key))
                            myfile.write('%s\tPacBio\tCDS\t%s\t%s\t.\t%s\t.\tgene_id "%s"; transcript_id "%s";\n'%(value[1], CDSexonSite[N],CDSexonSite[N+1],value[0],value[2],key))
                            N+=2
                    else:#位点不是成对的
                        #print len(CDSexonSite),CDSexonSite
                        print(key,'?2'*1000)

        else:#没有内含子，直接进行计算#完成
            #进行判断在CDS中有没有结果，有结果的进行计算
            if key in CDSsite.keys():#如果优结果
                #print key
                if value[0]=='+':
                    CDSstartSite=min(value[3])+CDSsite[key][0]-1
                    CDSstopSite=CDSstartSite+CDSsite[key][2]+3#为起始位置加上CDS长度
                    myfile.write('%s\t%s\tCDS\t%s\t%s\t.\t%s\t.\tgene_id "%s"; transcript_id "%s";\n'%(value[1],'PacBio',CDSstartSite,CDSstopSite,value[0],value[2],key))
                    #print key,value,CDSsite[key]
                if value[0]=='-':
                    #print key,value,CDSsite[key]
                    CDSstartSite=max(value[3])-CDSsite[key][0]+1
                    CDSstopSite=CDSstartSite-CDSsite[key][2]-3#终止位点为起始位点减去CDS长度
                    #print CDSstartSite,CDSstopSite
                    myfile.write('%s\t%s\tCDS\t%s\t%s\t.\t%s\t.\tgene_id "%s"; transcript_id "%s";\n'%(value[1],'PacBio',CDSstopSite,CDSstartSite,value[0],value[2],key))
    myfile.close()
start=time.time()
#FGid('max-CDS.fasta','all_tissue.collapsed.gff','tran-intron.gff')
FGid(args.cds,args.gtf,args.intron,args.out)
end=time.clock()
#print('Running time : %s '%(end-start))
