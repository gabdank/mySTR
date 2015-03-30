his = {}
bedFile = open("/media/gabdank/Disk3/mySTR/WS243.one.bed","r")
for l in bedFile:
    #print l
    arrL = l.strip().split()
    chromosome = arrL[0]
    start = int(arrL[1])
    end = int(arrL[2])
    size = end-start+1
    if not size in his:
        his[size]=0
    his[size]+=1

        
    #break
bedFile.close()
summa = 0
for size in sorted(his.keys()):
    summa += his[size]
    
    print str(size)+"\t"+str(100*(float(summa)/3182.0))