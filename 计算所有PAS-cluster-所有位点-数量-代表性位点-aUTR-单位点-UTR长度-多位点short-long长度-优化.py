#coding:utf-8
#2022.9.2
#luping
#该脚本用于计算APA的不同single short long 3UTR长度。以及aUTR的长度

def FGid(APAcount, countNum, LocationClass):
    # 1 第一步 先选取每个基因的主要表达cluster代表位点#cluster需要多个部分补充#目前APA是有两个位点的进行cluster，需要补充单个FLNC的结果#此外还有部分Isogene没有获得FLNC需要额外
    myfile = open('Allsite-MainSite2.txt', 'w')
    myfile.write('ID\tAll\tMian\n')
    myfile2 = open('Single-Short-long-3UTR2.txt', 'w')
    myfile2.write('ID\tType\tLengt\n')
    myfile3 = open('aUTR2.txt', 'w')
    myfile3.write('ID\tMaxSite\tMinSite\taUTR\tStrand\n')
    GeneTotal ={}#存放基因信息
    for line1 in open(APAcount):
        if 'PB' in line1:
            if 'OrginT' not in line1:
                FGid, CHR, Strand, APAnum, Site = line1.strip().split('\t')
                GeneTotal.setdefault(FGid,{})
                GeneTotal[FGid]['Chr'] = CHR
                GeneTotal[FGid]['Strand'] = Strand
                GeneTotal[FGid]['APAnum'] = int(APAnum)
                SiteRe = []  # 把位点信息还原成列表
                SiteReAddNum = []  # 位点信息加num统计
                SiteReAddNumClassT = []  # 位点信息加num统计的主要表达位点
                Site = Site.replace('[[','[').replace(']]',']').replace('\n','')
                # print('Site',Site)
                Location = Site
                Location2 = Location.split(']')
                # print(Location2)
                Location3 = []
                for i in Location2:
                    Cluster = i.replace('[', '').split(', ')
                    Cluster = list(filter(None, Cluster))
                    # print('Cluster',Cluster)
                    Clusterlist = []
                    ClusterlistNum = []  # 位点数量#从此处选出代表位点，进行UTR计算
                    for i in Cluster:  # 输出所有位点和数量文件加代表性位点和数量#输出代表性位点，单个UTR长度#输出多个代表性位点是否3UTR并计算aUTR和3UTR长度
                        # print(i)
                        Clusterlist.append(int(i))  #
                    Location3.append(Clusterlist)
                Location3 = list(filter(None,Location3))
                GeneTotal[FGid]['Location'] = Location3

    APAnumFile = {}  # 存放基因信息
    for line2 in open(countNum):
        PBid, PBsite, SiteNumber = line2.strip().split('\t')
        APAnumFile.setdefault(PBsite, {})
        APAnumFile[PBsite]['SiteNumber'] = int(SiteNumber)

    LocationClassFile = {}#
    for line3 in open(LocationClass):
        if 'PB' in line3:
         if '3UTR' in line3:
            Gene,ClassT,PASs,PAS,PasType,Strand,CDSrange,UTR5range,UTR3range,IntronRange,TranscriptRange = line3.strip().split('\t')
            LocationClassFile.setdefault(Gene,{})
            LocationClassFile[Gene]['ClassT'] = ClassT
            LocationClassFile[Gene]['PASs'] = PASs
            # LocationClassFile[Gene]['PAS'] = PAS
            # LocationClassFile[Gene]['PASType'] = PasType
            if 'PASType' not in LocationClassFile[Gene].keys():#需要3UTR类型的转录本，，将类型和位点添加
                LocationClassFile[Gene]['PASType'] = [PasType]
                LocationClassFile[Gene]['PAS'] = [PAS]
            else:
                PASTypeValue = LocationClassFile[Gene]['PASType']
                PASTypeValue.append(PasType)
                PASsValue = LocationClassFile[Gene]['PAS']
                PASsValue.append(PAS)
            LocationClassFile[Gene]['Strand'] = Strand
            LocationClassFile[Gene]['CDSrange'] = CDSrange
            LocationClassFile[Gene]['UTR5range'] = UTR5range
            LocationClassFile[Gene]['UTR3range'] = UTR3range
            LocationClassFile[Gene]['IntronRange'] = IntronRange
            LocationClassFile[Gene]['TranscriptRange'] = TranscriptRange

    #1 需要将位点与其数量对接，选出每簇的代表性位点
    for key,value in GeneTotal.items():
        # print(key,value['Location'])
        #遍历位点，选出支持位点最大
        ClusterName = []#位点名称#位点统计
        ClusterNumMax = []#最大位点统计
        for i in value['Location']:#遍历簇
            # print(key,value['Location'],i)
            EverCluster=[];EverClusterNum=[]#簇位点每类的名称和数量
            MaxEverCluster=[];MaxEverClusterNum=[]#簇位点每类的名称和数量,最大
            for j in i:#遍历簇里每个位点
                # print(key,i,j)
                #此处根据基因名，位点。索引第二个字典，收集位点数量
                GeneSite = '%s-%s' % (key, j)
                # print(GeneSite)
                EverCluster.append(GeneSite)
                EverClusterNum.append(APAnumFile[GeneSite]['SiteNumber'])
            # print(key,i,EverCluster,EverClusterNum,max(EverClusterNum),EverCluster[EverClusterNum.index(max(EverClusterNum))])
            MaxEverCluster.append((EverCluster,EverClusterNum))#所有位点记录
            MaxEverClusterNum.append((EverCluster[EverClusterNum.index(max(EverClusterNum))],str(max(EverClusterNum))))#最大位点记录
            # print(MaxEverCluster,MaxEverClusterNum)
            ClusterName.append(MaxEverCluster)
            ClusterNumMax.append(MaxEverClusterNum)#每簇加入簇列表
        # print(key,ClusterName,'LS'*10,ClusterNumMax)
        myfile.write(key+'\t'+str(ClusterName)+'\t'+str(ClusterNumMax)+'\n')
    #2 计算aUTR 和single，short，long 3UTR的长度
    for key3,value3 in LocationClassFile.items():
        if value3['CDSrange'] != '[]':#转录本必须有CDS#
            # print(key3,value3['PAS'],value3['Strand'],value3['CDSrange'])
            CDS = value3['CDSrange'].replace('[','').replace(']','').replace("'","").replace('-', ', ').split(', ')
            # print(CDS)
            CDS =[int(i) for i in CDS]
            CDSMinMax =[];CDSMinMax.append(min(CDS));CDSMinMax.append(max(CDS))
            PAS = value3['PAS']
            PAS = [int(j) for j in PAS]
            # print(key3,PAS)
            #print(key3+'\t'+str(max(PAS))+'\t'+str(min(PAS))+'\t' +str(max(PAS)-min(PAS))+'\t'+value3['Strand']+'\n')
            if len(value3['PAS'])>1:#如果3UTR位点有两个及以上，可计算aUTR，short，long
                myfile3.write(key3+'\t'+str(max(PAS))+'\t'+str(min(PAS))+'\t' +str(max(PAS)-min(PAS))+'\t'+value3['Strand']+'\n')
                if value3['Strand'] =='+':
                    ShortUTR = min(PAS) - max(CDS)
                    LongUTR = max(PAS) - max(CDS)
                    # print(key3+'\t'+'ShortLest'+'\t'+str(ShortUTR)+'\n')
                    # print(key3+'\t'+'Longest'+'\t'+str(LongUTR)+'\n')
                    myfile2.write(key3+'\t'+'ShortLest'+'\t'+str(ShortUTR)+'\n')
                    myfile2.write(key3+'\t'+'Longest'+'\t'+str(LongUTR)+'\n')
                else:
                    ShortUTR = min(CDS) - max(PAS)
                    LongUTR = min(CDS) - min(PAS)
                    myfile2.write(key3 + '\t' + 'ShortLest' + '\t' + str(ShortUTR) + '\n')
                    myfile2.write(key3 + '\t' + 'Longest' + '\t' + str(LongUTR) + '\n')
            if len(value3['PAS'])== 1:  #如果3UTR位点只有一个，计算为single 3UTR，不计算aUTR
                if value3['Strand'] =='+':
                    SingleUTR = min(PAS)-max(CDS)
                    myfile2.write(key3 + '\t' + 'Single' + '\t' + str(SingleUTR) + '\n')
                if value3['Strand'] == '-':
                    SingleUTR = min(CDS)-min(PAS)
                    myfile2.write(key3 + '\t' + 'Single' + '\t' + str(SingleUTR) + '\n')
    myfile.close()
    myfile2.close()
    myfile3.close()

FGid('all_tissue-tama.tracking-PAS-number-count.txt-polyA-site-count.txt', 'all_tissue-tama.tracking-PAS-number-count.txt',
     'PAS-location-count.txt')
