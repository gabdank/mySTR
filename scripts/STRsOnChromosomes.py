import numpy as np
import matplotlib.pyplot as plt
chrLengthDict = {"I":15072423, "II":15279345,"III":13783700,"IV":17493793,"V":20924149,"X":17718866,"M":13794, "chrMtDNA":13794}
coordinated={'I':(3858000,11040000),'II':(4879000,12020000),'III':(3722000,10340000),'IV':(3896000,12970000),'V':(5897000,16550000),'X':(6137000,12480000)}
his = {}
binSize = 50000
chromosomeName = 'X'
bedFile = open("/media/gabdank/Disk3/mySTR/WS243.one.bed","r")
for l in bedFile:
    #print l
    arrL = l.strip().split()
    chromosome = arrL[0]
    start = int(arrL[1])
    end = int(arrL[2])
    size = end-start+1
    center = (size/2)+start
    binnedCoordinate = center/binSize    
    
    if chromosome == chromosomeName:
        if not binnedCoordinate in his:
            his[binnedCoordinate]=0
        his[binnedCoordinate]+=1
        
bedFile.close()

x_ticks = np.arange(0,chrLengthDict[chromosomeName],5000000)
x_v = []
y_v = []
for size in sorted(his.keys()):
    x_v.append(size*binSize)
    y_v.append(his[size])

plt.bar(x_v,y_v,100000)
plt.axvline(x=coordinated[chromosomeName][0], ymin=0, ymax = 30, linewidth=2, color='red')
plt.axvline(x=coordinated[chromosomeName][1], ymin=0, ymax = 30, linewidth=2, color='red')
plt.grid()
plt.xticks(x_ticks)
plt.title('STRs distribution across chr'+chromosomeName+" bin size 100K")
plt.ylabel('Number of occurrences')
plt.xlabel('Position across the chromosome')
plt.xlim([0,chrLengthDict[chromosomeName]])
plt.show()
