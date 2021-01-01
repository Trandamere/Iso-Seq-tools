#coding:utf-8
def FGid(gff):
    myfile1=open('PAS-Up-60-25bp.fasta','w')
    myfile2=open('PAS-Up-25-15bp.fasta','w')
    myfile3=open('PAS-Up-15-0bp.fasta','w')
    myfile4=open('PAS-Down-0-15bp.fasta','w')
    myfile5=open('PAS-Down-15-60bp.fasta','w')
    myfile6=open('PAS-Up-25-0bp.fasta','w')
    title=[]
    fasta=[]
    for line in open(gff):
        if '>' in line:
            title.append(line)
        else:
            fasta.append(line)
    while len(title)>0:
        titlepop=title.pop()
        fastapop=fasta.pop()
        #print titlepop
        #print fastapop[39:79]#Up60-U25#35长度向后延伸5bp为了6聚体统计
        myfile1.write(titlepop)
        myfile1.write('%s\n'%fastapop[39:79])
        #print titlepop
        #print fastapop[74:89]#Up25-U15#10长度向后延伸5bp为了6聚体统计
        myfile2.write(titlepop)
        myfile2.write('%s\n'%fastapop[74:89])        
        #print titlepop
        #print fastapop[84:104]#Up15-0#长度向后延伸5bp为了6聚体统计
        myfile3.write(titlepop)
        myfile3.write('%s\n'%fastapop[84:104])
        myfile6.write(titlepop)
        myfile6.write('%s\n'%fastapop[74:104])        
        #print titlepop
        #print fastapop[100:120]#下游0-15#长度向后延伸5bp为了6聚体统计
        myfile4.write(titlepop)
        myfile4.write('%s\n'%fastapop[100:120])
        #print titlepop
        #print fastapop[115:165]#下游15-60#长度向后延伸5bp为了6聚体统计
        myfile5.write(titlepop)
        myfile5.write('%s\n'%fastapop[115:165])
    myfile1.close()
    myfile2.close()
    myfile3.close()
    myfile4.close()
    myfile5.close()
    myfile6.close()
#FGid('join-PAS-updown201bp.fasta')
FGid('join-PAS-updown201bp-all.fasta')
