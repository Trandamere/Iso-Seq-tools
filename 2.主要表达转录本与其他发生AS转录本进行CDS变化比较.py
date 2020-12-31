#coding:utf-8
#2020.7.31
#根据主要表达转录本，进行其他亚型比较
def CDSchange(GeneT,ASgtf,PH1,Allgtf):
    myfile=open('ASgene-CDS-change-count.txt','w')
    myfile.write('Gene\tRefTranscript\tAStranscript\tCDSType\n')
    for line1 in open(GeneT):
        if 'FG' in line1:
            ASgene=line1.split('\t')[0]
            Tone=line1.strip().split('\t')[-1]
            #print ASgene,T
            AST=[]#收集发生AS的转录本
            for line2 in open(ASgtf):
                if '%sT'%ASgene in line2:
                    AS1event=line2.split('transcript1_id "')[1].split('"')[0]
                    AS2event=line2.split('transcript2_id "')[1].split('"')[0]
                    if AS1event!=Tone:#不添加主要转录本
                        AST.append(AS1event)
                    if AS2event!=Tone:#不添加主要转录本#下后续去除外显子结构相同的转录本，来排除相同剪切类型转录本
                        AST.append(AS2event)
            #先收集比较AS转录本exon，去除
            Tsj=[];ASTsj=[]#收集主要转录本和AS转录本的SJ位点，进行比较
            for line3 in open(PH1):
                if '%s"'%Tone in line3:#收集主要转录本的剪切位点
                    if 'exon\t' in line3:
                        #print ASgene,T,line3
                        Tsj.append(int(line3.split('\t')[3]))
                        Tsj.append(int(line3.split('\t')[4]))
            Tsj.sort()#代表转录本位点排序
            del Tsj[0]
            del Tsj[-1]#只保留主要转录本内部剪切位点
            
            ASTA=[]
            for AA in AST:
                if AA not in ASTA:
                    ASTA.append(AA)
            ASTA.sort()
            #print len(ASTA)
            #print ASgene,Tone,ASTA
            Rmove=[]
            for I in ASTA:#收集AS转录本的剪切位点###为什么for循环会少最后一个
                #print I,'I'
                Tallsj=[]#单个转录本所有外显子位置
                for line4 in open(PH1):
                    if '%s"'%I in line4:
                        if 'exon\t' in line4:
                            #print ASgene,T,I,line4
                            Tallsj.append(int(line4.split('\t')[3]))
                            Tallsj.append(int(line4.split('\t')[4]))
                            
                Tallsj.sort()#AS转录本位点排序
                del Tallsj[0]
                del Tallsj[-1]
                #print ASgene,Tone,Tsj,I,Tallsj,'Same'*10
                if Tsj==Tallsj:
                    #i是不符合的转录本
                    #print ASgene,T,AST
                    #print 'Gene modle',T,'include',i
                    Rmove.append(I)#去掉SJ位点相同类型的转录本
            #print 'Adjest',ASgene,T,AST
            for res in Rmove:#单独删除#直接在列表内删除会引发bug
                ASTA.remove(res)
            #print 'Afer move',Rmove,ASTA
            #part2 接下来进行主要转录本和AS转录本比较，进行CDS分类
            #T为参考转录本，AST为过滤后的AS转录本#增加CDS延长类型统计
            Ta=Tone
            for Tb in ASTA:
                strand=[];mRNATa=[];mRNATb=[];CDSTa=[];CDSTb=[]#收集转录本起始终止位点和CDS起始终止位点
                for line5 in open(Allgtf):
                    if 'CDS\t' in line5:
                        if '%s"'%Ta in line5:
                            strand.append(line5.split('\t')[6])
                            CDSTa.append(int(line5.split('\t')[3]))
                            CDSTa.append(int(line5.split('\t')[4]))
                        if '%s"'%Tb in line5:
                            CDSTb.append(int(line5.split('\t')[3]))
                            CDSTb.append(int(line5.split('\t')[4]))
                    if 'transcript\t' in line5:
                        if '%s"'%Ta in line5:
                            mRNATa.append(int(line5.split('\t')[3]))
                            mRNATa.append(int(line5.split('\t')[4]))
                        if '%s"'%Tb in line5:
                            mRNATb.append(int(line5.split('\t')[3]))
                            mRNATb.append(int(line5.split('\t')[4]))
                CDSTa=sorted(CDSTa);CDSTb=sorted(CDSTb)#排序
                #区分是否是ATI和APA不分链方向
                #print Ta,Tb,mRNATa,mRNATb,CDSTa,CDSTb
                if len(CDSTb) and len(mRNATb) and len(CDSTa) and len(mRNATa)>0:
                    if min(mRNATa)<=min(CDSTb) and min(mRNATb)<=min(CDSTa) and max(mRNATa)>=max(CDSTb) and max(mRNATb)>=max(CDSTa):#保证不受ATI和APA干扰
                        if strand[0]=='+':
                            UTR5=min(CDSTb)-min(CDSTa)
                            UTR3=max(CDSTb)-max(CDSTa)
                            if UTR5>0:#AS转录本CDS起始靠后
                                if UTR3>0:#AS转录本CDS终止靠后
                                    print ASgene,Ta,Tb,"CDS moved"
                                    myfile.write('%s\t%s\t%s\t%s\n'%(ASgene,Ta,Tb,"CDS moved"))
                                elif UTR3==0:#AS转录本CDS终止相同
                                    print ASgene,Ta,Tb,"Only 5'CDS short"
                                    myfile.write('%s\t%s\t%s\t%s\n'%(ASgene,Ta,Tb,"Only 5'CDS short"))
                                elif UTR3<0:#AS转录本CDS终止靠前
                                    print ASgene,Ta,Tb,"5'CDS and 3'CDS short"
                                    myfile.write('%s\t%s\t%s\t%s\n'%(ASgene,Ta,Tb,"5'CDS and 3'CDS short"))
                                else:
                                    print ASgene,Ta,Tb,'?'*100
                            if UTR5==0:#AS转录本CDS起始相同
                                if UTR3>0:#AS转录本CDS终止靠后
                                    print ASgene,Ta,Tb,"Only 3'CDS long"
                                    myfile.write('%s\t%s\t%s\t%s\n'%(ASgene,Ta,Tb,"Only 3'CDS long"))
                                elif UTR3==0:#AS转转录本终止位点相同
                                    if CDSTa==CDSTb:#若CDS区域完全一样，说明AS发生在UTR
                                        print ASgene,Ta,Tb,"CDS maintain"
                                        myfile.write('%s\t%s\t%s\t%s\n'%(ASgene,Ta,Tb,"CDS maintain"))
                                    else:
                                        print ASgene,Ta,Tb,"3n in CDS"#AS区域是3n
                                        myfile.write('%s\t%s\t%s\t%s\n'%(ASgene,Ta,Tb,"3n in CDS"))
                                elif UTR3<0:##AS转录本CDS终止靠前
                                    print ASgene,Ta,Tb,"Only 3'CDS short"#AS区只有3CDS缩短
                                    myfile.write('%s\t%s\t%s\t%s\n'%(ASgene,Ta,Tb,"Only 3'CDS short"))
                                else:
                                    print ASgene,Ta,Tb,'?'*100
                            if UTR5<0:#AS转录本CDS起始变长
                                if UTR3>0:#AS转录本CDS终止靠后
                                    print ASgene,Ta,Tb,"5'CDS and 3'CDS long"
                                    myfile.write('%s\t%s\t%s\t%s\n'%(ASgene,Ta,Tb,"5'CDS and 3'CDS long"))
                                elif UTR3==0:#AS转录本CDS终止相同
                                    print ASgene,Ta,Tb,"Only 5'CDS long"
                                    myfile.write('%s\t%s\t%s\t%s\n'%(ASgene,Ta,Tb,"Only 5'CDS long"))
                                elif UTR3<0:#AS转录本CDS终止靠前
                                    print ASgene,Ta,Tb,"CDS moved"
                                    myfile.write('%s\t%s\t%s\t%s\n'%(ASgene,Ta,Tb,"CDS moved"))
                                else:
                                    print ASgene,Ta,Tb,'?'*100
                        if strand[0]=='-':
                            UTR3=min(CDSTa)-min(CDSTb)
                            UTR5=max(CDSTa)-max(CDSTb)
                            if UTR5>0:#AS转录本CDS起始靠后
                                if UTR3>0:#AS转录本CDS终止靠后
                                    print ASgene,Ta,Tb,"CDS moved"
                                    myfile.write('%s\t%s\t%s\t%s\n'%(ASgene,Ta,Tb,"CDS moved"))
                                elif UTR3==0:#AS转录本CDS终止相同
                                    print ASgene,Ta,Tb,"Only 5'CDS short"
                                    myfile.write('%s\t%s\t%s\t%s\n'%(ASgene,Ta,Tb,"Only 5'CDS short"))
                                elif UTR3<0:#AS转录本CDS终止靠前
                                    print ASgene,Ta,Tb,"5'CDS and 3'CDS short"
                                    myfile.write('%s\t%s\t%s\t%s\n'%(ASgene,Ta,Tb,"5'CDS and 3'CDS short"))
                                else:
                                    print ASgene,Ta,Tb,'?'*100
                            if UTR5==0:#AS转录本CDS起始相同
                                if UTR3>0:#AS转录本CDS终止靠后
                                    print ASgene,Ta,Tb,"Only 3'CDS long"
                                    myfile.write('%s\t%s\t%s\t%s\n'%(ASgene,Ta,Tb,"Only 3'CDS long"))
                                elif UTR3==0:#AS转转录本终止位点相同
                                    if CDSTa==CDSTb:#若CDS区域完全一样，说明AS发生在UTR
                                        print ASgene,Ta,Tb,"CDS maintain"
                                        myfile.write('%s\t%s\t%s\t%s\n'%(ASgene,Ta,Tb,"CDS maintain"))
                                    else:
                                        print ASgene,Ta,Tb,"3n in CDS"#AS区域是3n
                                        myfile.write('%s\t%s\t%s\t%s\n'%(ASgene,Ta,Tb,"3n in CDS"))
                                elif UTR3<0:##AS转录本CDS终止靠前
                                    print ASgene,Ta,Tb,"Only 3'CDS short"#AS区只有3CDS缩短
                                    myfile.write('%s\t%s\t%s\t%s\n'%(ASgene,Ta,Tb,"Only 3'CDS short"))
                                else:
                                    print ASgene,Ta,Tb,'?'*100
                            if UTR5<0:#AS转录本CDS起始变长
                                if UTR3>0:#AS转录本CDS终止靠后
                                    print ASgene,Ta,Tb,"5'CDS and 3'CDS long"
                                    myfile.write('%s\t%s\t%s\t%s\n'%(ASgene,Ta,Tb,"5'CDS and 3'CDS long"))
                                elif UTR3==0:#AS转录本CDS终止相同
                                    print ASgene,Ta,Tb,"Only 5'CDS long"
                                    myfile.write('%s\t%s\t%s\t%s\n'%(ASgene,Ta,Tb,"Only 5'CDS long"))
                                elif UTR3<0:#AS转录本CDS终止靠前
                                    print ASgene,Ta,Tb,"CDS moved"
                                    myfile.write('%s\t%s\t%s\t%s\n'%(ASgene,Ta,Tb,"CDS moved"))
                                else:
                                    print ASgene,Ta,Tb,'?'*100
                    else:#起始和终止处在CDS区域
                        print ASgene,Ta,Tb,'ATI or APA'
                        myfile.write('%s\t%s\t%s\t%s\n'%(ASgene,Ta,Tb,'ATI or APA'))
                else:
                    print ASgene,Ta,Tb,'CDS vanished'
                    myfile.write('%s\t%s\t%s\t%s\n'%(ASgene,Ta,Tb,'CDS vanished'))
    myfile.close()
  
#CDSchange('Gene-T-for-AS.txt','All-landscape.gtf','Fusion-PH1-fusion-version2.3.gtf.gff3.gtf','PH1-fusion-version3.3.gtf')
#CDSchange('Gene-T-for-AS2.txt','All-landscape.gtf','Fusion-PH1-fusion-version2.3.gtf.gff3.gtf','PH1-fusion-version3.3.gtf')
CDSchange('Gene-T-for-AS.txt','delMt-All-landscape.gtf','delMt-Fusion-PH1-fusion-version2.3.gtf.gff3.gtf','PH1-fusion-version4.1.gtf')
