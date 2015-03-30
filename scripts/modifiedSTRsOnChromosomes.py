import numpy as np
import matplotlib.pyplot as plt
chrLengthDict = {"I":15072423, "II":15279345,"III":13783700,"IV":17493793,"V":20924149,"X":17718866,"M":13794, "chrMtDNA":13794}
coordinated={'I':(3858000,11040000),'II':(4879000,12020000),'III':(3722000,10340000),'IV':(3896000,12970000),'V':(5897000,16550000),'X':(6137000,12480000)}

def interesting(d,strain):   
    if d[(strain,'lower')]!='0':
        return True
    else:
        refLength = d[('reference','length')]
        strainLength = d[(strain,'length')]
        strainSTD = d[(strain,'STD')]
        if strainLength!=0 and (refLength<strainLength-(2*strainSTD) or refLength>strainLength+(2*strainSTD)):
            return True
    return False
            
def gatherInformation(l):
    d = {}
    arr = l.strip().split('\t')
    referenceLength = int(arr[1].split(']')[1])    
    d[('reference','length')]=referenceLength
    d[('MY1','length')]=float(arr[2])
    d[('MY1','STD')]=float(arr[3])
    d[('MY1','lower')]=arr[4]
    d[('MY2','length')]=float(arr[5])
    d[('MY2','STD')]=float(arr[6])
    d[('MY2','lower')]=arr[7]
    d[('AB1','length')]=float(arr[8])
    d[('AB1','STD')]=float(arr[9])
    d[('AB1','lower')]=arr[10]
    d[('AB3','length')]=float(arr[11])
    d[('AB3','STD')]=float(arr[12])
    d[('AB3','lower')]=arr[13]    
    d[('MY6','length')]=float(arr[14])
    d[('MY6','STD')]=float(arr[15])
    d[('MY6','lower')]=arr[16]
    d[('MY14','length')]=float(arr[17])
    d[('MY14','STD')]=float(arr[18])
    d[('MY14','lower')]=arr[19]
    d[('MY16','length')]=float(arr[20])
    d[('MY16','STD')]=float(arr[21])
    d[('MY16','lower')]=arr[22]
    return d

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
    
mutatedSTRs = {}
mone = 0
mergedFile = open("outputFile","r")
mergedFile.readline()

chromoso = chromosomeName
strainNames = ['MY1','MY2','AB1','AB3','MY6', 'MY14', 'MY16']
for strain in strainNames:
    mutatedSTRs[strain]=[]

for l in mergedFile:    
    arr = l.split("\t")
    chromo = arr[0].split(':')[0]
    position = arr[0].split(':')[1]
    d = gatherInformation(l)
    for strain in strainNames:   
        if interesting(d,strain) and chromo==chromoso:            
            mutatedSTRs[strain].append(position)
        #print chromo+":"+str(position)+"\t"+str(d['reference','length'])+"\t"+ str(d[(strain,'length')])+"\t"+str(d[(strain,'lower')])     
mergedFile.close()
#print mutatedSTRs['MY1']
    

plt.bar(x_v,y_v, 100000)
plt.axvline(x=coordinated[chromosomeName][0], ymin=0, ymax = 40, linewidth=2, color='red')
plt.axvline(x=coordinated[chromosomeName][1], ymin=0, ymax = 40, linewidth=2, color='red')

mutVals = {}
appendableValue = 0
for mut in mutatedSTRs.keys():
    mutVals[mut]=[]
    appendableValue -= 1
    for v in mutatedSTRs[mut]:
        mutVals[mut].append(appendableValue)    
    plt.scatter(mutatedSTRs[mut],mutVals[mut])
#for v in mutatedSTRs['MY6']:
   
#print len(mutVals)
#print len(mutatedSTRs['MY6'])
   
#for v in mutatedSTRs['MY16']:
#    plt.axvline(x=v, ymin=0, ymax = 40, linewidth=1, color='yellow')
    


plt.grid()
plt.xticks(x_ticks)
plt.title('STRs distribution across chr'+chromosomeName)
plt.ylabel('Number of occurrences')
plt.xlabel('Position across the chromosome')
plt.xlim([0,chrLengthDict[chromosomeName]])
plt.show()
