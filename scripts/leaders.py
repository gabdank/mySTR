def readFile(fileName, numOfLeaders):
    leadersDict = {}
    fileRead = open(fileName,"r")
    fileRead.readline()
    mone = 0
    for l in fileRead:
        mone += 1
        arr = l.split()
        chr = arr[0].split(':')[0]
        pos = arr[0].split(':')[1]
        key = (chr,pos)
        leadersDict[key]=1
        
        if mone>=numOfLeaders:
            break
    
    fileRead.close()
    return leadersDict

my1D = readFile("my1.ordered",40)
my2D = readFile("my2.ordered",40)
my6D = readFile("my6.ordered",40)
ab1D = readFile("ab1.ordered",40)

for x in my1D:
    z = 0    
    if x in my2D: 
        z += 1
    if x in my6D:
        z += 1
    if x in ab1D: 
        z += 1
    if z>1:
        print x
            

