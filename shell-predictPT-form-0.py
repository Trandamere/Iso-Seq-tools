#coding:utf-8

#2018.12.20
#luping
#从原始矫正的fasta序列中选取多顺反子
#输入文件为原始序列生成的gff3
import os

#函数将gff3生成gtf
def gff3togtf(x):
    myfile=open(r'tofu-%s.gtf'%x,'w')#
    for line in open(r'%s'%x):
        if 'mRNA\t' in line:
            a=line.split('ID')[0]#注释前半部分
            tranid=(line.split('=')[1]).split('.mrna')[0]
            myfile.write('%sgene_id "%s"; transcript_id "%s";\n'%(a,tranid,tranid))
        if 'exon\t' in line:
            a=line.split('ID')[0]#注释前半部分
            tranid=(line.split('=')[1]).split('.mrna')[0]
            myfile.write('%sgene_id "%s"; transcript_id "%s";\n'%(a,tranid,tranid))
    myfile.close()
    
def findPT(x):
    tranidlist=[]
    for line in open(r'ORF-%s'%x):
        if 'mRNA\t' in line:
            tranidlist.append((line.split('.p')[0]).split('=')[1])
    tranid=list(set(tranidlist))
    myfile=open(r'%s-polycistronic-transcript.txt'%(x.split('-')[0]),'w')
    myfile2=open(r'%s-polycistronic-transcript+start-stop.txt'%(x.split('-')[0]),'w')
    while len(tranid)>0:
        print 'remain tran %s'%len(tranid)
        SS=[]#存放ORF的区域
        pop=tranid.pop()
        tran=pop#转录本号备份
        for line in open(r'ORF-%s'%x):
            if 'mRNA\t' in line:
                if '%s.p'%pop in line:
                    A=[int(line.split('\t')[3]),int(line.split('\t')[4])]#提取基因号所属的ORF区域
                    SS.append(A)#将所有的ORF区域添加到SS中
        SS2=SS#复制SS生成相同的ORF编码区域集合
        SSall=[]#将ORF位置进行两两组合
        for i in SS:
            for j in SS2:
                if i==j:
                    pass
                else:
                    SSall.append([i,j])#将基因下属的ORF进行亮亮组合
        SSover=[]#存放有重叠的ORF组合
        for i in SSall:
            for j in range(i[0][0],i[0][1]):
                if j in range(i[1][0],i[1][1]):
                    SSover.append(i)#将有区域重叠的mRNA组合提取，判断在两个ORF中是否有公共领域
        SSover2=[]#将无重复orf结果去重复
        for i in SSover:
            if i in SSover2:
                pass
            else:
                SSover2.append(i)
        final=[]#将有重复的ORF从全部集合中取出
        for i in SSover2:
            if i in SSall:
                SSall.remove(i)#将ORF组合逐个去除
        SSall2=[]
        for i in SSall:
            if i in SSall2:
                pass
            else:
                SSall2.append(i)
        if len(SSall)>0:#若最后的ORF组合中还有ORF组合剩余，就是有没有重叠ORF组合
            myfile.write('%s\n'%(tran))#记录有两个不重叠ORF的转录本
            myfile2.write('%s\t%s\n'%(tran,SSall2))#记录有两个不同重叠ORF转录本和编码区域

def PTgtf(x):
    PBtran=[]
    for line1 in open(r'%s-polycistronic-transcript.txt'%(x.split('-')[0])):
        PBtran.append(line1.strip())
    PBtran=list(set(PBtran))
    print '%s\tPBtran\t%s'%((x.split('-')[0]),len(PBtran))
    myfile1=open(r'%s-polycistronic-transcript.gtf'%(x.split('-')[0]),'w')
    for line2 in open(r'ORF-%s'%x):
        if 'mRNA\t' in line2:
            if (line2.split('.p')[0]).split('=')[1] in PBtran:
                myfile1.write(line2)
        if 'exon\t' in line2:
            if (line2.split('.p')[0]).split('=')[1] in PBtran:
                myfile1.write(line2)
    myfile1.close()

def CompareToFG(x,y):
    os.system(r'/home/luping/gffcompare-0.10.2/gffcompare -r %s-polycistronic-transcript.gtf -o ./%s-%s-polycistronic-transcript-rev %s'%(x.split('-')[0],x.split('-')[0],y.split('.')[0],y))
    os.system(r'/home/luping/gffcompare-0.10.2/gffcompare -r %s -o ./%s-%s-polycistronic-transcript %s-polycistronic-transcript.gtf'%(y,y.split('.')[0],x.split('-')[0],x.split('-')[0]))

def PTfromGffcompare(x,y):
    same=['=','c','k','j','e','o']
    tranid=[]
    for line in open(r'%s-polycistronic-transcript.txt'%(x.split('-')[0])):
        tranid.append(line.strip())
    tranid2=list(set(tranid));tranid3=list(set(tranid2))#备份正反向消耗
    myfile3=open(r'%s-PT-1-n.txt'%(y.split('.')[0]),'w')
    while len(tranid2)>0:
        print len(tranid2)
        FGPB=[]
        pop=tranid2.pop()
        for line in open(r'%s-%s-polycistronic-transcript-rev.tracking'%(x.split('-')[0],y.split('.')[0])):
            if '%s.p'%pop in line:
                if line.split('\t')[3] in same:
                    a=[(line.split('\t')[2]).split('.p')[0],((line.split('\t')[4]).split('|')[0]).split(':')[1]]#提取FG与PB
                    b='%s\t%s'%(a[0],a[1])#修改PB-FG输出模式
                    FGPB.append(b)
        FGPB2=list(set(FGPB))
        if len(FGPB)<2:
            pass
        elif FGPB2[0].split('\t')[0]==FGPB2[1].split('\t')[0]:#将转录本下的cds比对结果集合，选出有2个以上的基因，若不同FG属于同一个PB，为目标
            for i in FGPB2:
                myfile3.write('%s\n'%i)

    while len (tranid3)>0:
        print len(tranid3)
        FGPB=[]
        pop=tranid3.pop()
        for line in open(r'%s-%s-polycistronic-transcript.tracking'%(y.split('.')[0],x.split('-')[0])):
            if '%s.p'%pop in line:
                if line.split('\t')[3] in same:
                    a=[(((line.split('\t')[4]).split('|')[0]).split(':')[1]).split('.p')[0],(line.split('\t')[2]).split('|')[0]]#提取FG与PB
                    b='%s\t%s'%(a[0],a[1])
                    FGPB.append(b)
        FGPB2=list(set(FGPB))
        if len(FGPB2)<2:
            pass
        elif FGPB2[0].split('\t')[0]==FGPB2[1].split('\t')[0]:
            for i in FGPB2:
                myfile3.write('%s\n'%i)
    myfile3.close()
    PTtran=[];PTgene=[]
    for line in open(r'%s-PT-1-n.txt'%(y.split('.')[0])):
        PTtran.append(line.split('\t')[0])
        PTgene.append(line.split('.')[1])
    PTtran=list(set(PTtran));PTgene=list(set(PTgene))
    print 'FGmapping tran\t%s'%(len(PTtran))
    print 'FGmapping gene\t%s'%(len(PTgene))

def PBinfasta(x,y,z,a,b):

    PTtran=[]
    for line in open(r'%s-PT-1-n.txt'%(x.split('.')[0])):
        if line.split('\t')[0] in PTtran:
            pass
        else:
            PTtran.append(line.split('\t')[0])
    for line in open(r'%s-PT-1-n.txt'%(y.split('.')[0])):
        if line.split('\t')[0] in PTtran:
            pass
        else:
            PTtran.append(line.split('\t')[0])
    PTreads=[]
    for line in open(r'%s'%z):
        if line.split('\t')[0] in PTtran:
            PB=((line.strip()).replace(',','\t')).split('\t')[1:]
            for i in PB:
                PTreads.append(i)
    #接下来从所有文件中去除顺反子序列

    myfile4=open(r'PT-%s'%a,'w')
    myfile5=open(r'filter-PT-%s'%a,'w')
    title=[]
    fasta=[]
    for line in open(r'%s'%a):
        if '>' in line:
            title.append(line)
        else:
            fasta.append(line)
    while len(title)>0:
        print len(title)
        titlepop=title.pop()
        fastapop=fasta.pop()
        if (titlepop.split(' ')[0]).split('>')[1] in PTreads:
            myfile4.write(titlepop)
            myfile4.write(fastapop)
        else:
            myfile5.write(titlepop)
            myfile5.write(fastapop)
    myfile4.close()
    myfile5.close()

    myfile6=open(r'filter-PT-%s'%b,'w')
    myfile7=open(r'PT-%s'%b,'w')
    for line in open(r'%s'%b):
        if line.split('\t')[0] in PTreads:
            myfile7.write(line)
        else:
            myfile6.write(line)
    myfile6.close()
    myfile7.close()

#2018.12.21改变为使用各基因的nofusion版本进行多顺反子预测，将符合标准的PB选出，进一步选出相关的原始序列
 
def PTpredictFromPB(x,y,z,a,b,c):
    #1 提取所有序列碱基进行ORF预测
    os.system('gffread -w gene-tofu-%s.gtf.fasta -g ph1.fasta %s'%(x,x))
    #2 进行ORF预测
    os.system('TransDecoder.LongOrfs -S -t gene-tofu-%s.gtf.fasta'%x)
    #3 将预测的CDS 进行基因组比对
    os.system('/disk/luping/tools/gmap-2017-11-15/bin/gmap -D /disk/luping/bam/iso-bam/pb-cluster-hqlq-total/two-condition/ref-FG -d PH1 -f 2 -n 1 -t 15 --no-chimeras --min-intronlength 20 --max-intronlength-middle 4000 --max-intronlength-ends 4000 -z sense_force ./gene-tofu-%s.gtf.fasta.transdecoder_dir/longest_orfs.cds >ORF-%s'%(x,x))
    #4 进行含有无重叠ORF的基因
    findPT(x)
    #6 从原始注释文件中提取顺反子候选注释
    PTgtf(x)
    #7 将多顺反子reads与FGSG,FGRAM分别比对
    CompareToFG(x,y)
    CompareToFG(x,z)
    #8 寻找CDS匹配2个FG的转录本
    PTfromGffcompare(x,y)
    PTfromGffcompare(x,z)
    #9 根据PB转录本号，选取原始reads#并从原始reads的fasta中去除，
    PBinfasta(y,z,a,b,c)
#PTpredictFromPB('ISO-nofusion.collapsed.gff','gff-FGSG.gff3','gff-FGRAMPH1.gff3','ISO-nofusion.collapsed.group.txt','ISO.fasta','ISO.sort.sam')
#PTpredictFromPB('Hy-ISO-nofusion.collapsed.gff','gff-FGSG.gff3','gff-FGRAMPH1.gff3','Hy-ISO-nofusion.collapsed.group.txt','Hy.fasta','Hy-ISO.sort.sam')
#PTpredictFromPB('Sex-ISO-fusion.collapsed.gff','gff-FGSG.gff3','gff-FGRAMPH1.gff3','Sex-ISO-nofusion.collapsed.group.txt','Sex.fasta','Sex-ISO.sort.sam')
