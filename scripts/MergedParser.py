from operator import itemgetter
#from TableCreator import processRecord

def parseChromosomeIndexToName(index):
    return 'chrI'

def processReferenceLine(line):
    arr = line.strip().split("$")
    chromosomeName = ""
    position = -1
    unit = ""
    repLength = -1
    for rec in arr:
        internalArr = rec.split(":")
        if internalArr[0]=='chromosome':
            chromosomeName = parseChromosomeIndexToName(int(internalArr[1]))
        if internalArr[0]=='position':
            position = int(internalArr[1])
        if internalArr[0]=='repetitiveUnit':
            unit = internalArr[1]
        if internalArr[0]=='repetitiveLength':
            repLength = int(internalArr[1])
    keyTuple = (chromosomeName,position)
    valueTuple = (unit, repLength/len(unit))
    return (keyTuple,valueTuple)

def processAlignedLine(line):
    arr = line.strip().split("$")
    for rec in arr:
        internalArr = rec.split(":")
        if internalArr[0]=='leftSequence':
            leftFlank = internalArr[1]
        if internalArr[0]=='rightSequence':
            rightFlank = internalArr[1]
        if internalArr[0]=='repetitiveSequence':
            repetitive = internalArr[1]
    return (len(leftFlank),len(repetitive),len(rightFlank))

def getRange(listOfVals):
    if len(listOfVals)>0:
        minimum = min(listOfVals)
        maximum = max(listOfVals)
        return (minimum,maximum)
    else:
        return ()

def getLeadingValue(listOfVals):
    valDic = {}
    for v in listOfVals:
        if int(v) in valDic:
            valDic[int(v)] += 1
        else:
            valDic[int(v)] = 1
    max = 0
    for v in valDic:
        if valDic[v]>max:
            max = valDic[v]

    maxCounter = 0
    location = None
    for v in valDic:
        if max - valDic[v] < 2:
            maxCounter += 1
            location = v
    return location

def readInMergedFile(fileName):
    locationsList =[]
    dict = {}
    mergedFile = open(fileName,"r")
    mone = 0

    for line in mergedFile:

        if line[0:1]=='r':
            mone += 1
            t = processReferenceLine(line)
            #print t
            currentKey = t[0]
            locationsList.append(currentKey)
            currentRef = t[1]
            dict[(currentKey,"reference")]=currentRef
            dict[(currentKey,"lower")]=[]
            dict[(currentKey,"exact")]=[]
        if line[0:1]=='[':
            #print currentKey
            #print currentRef
            # potentially can do some refinment on the alignment techniques.. not sure it is necessary

            leftRepRight = processAlignedLine(line[1:-2])
            if (leftRepRight[0]<5 or leftRepRight[2]<5): #lowerBound
                dict[(currentKey,"lower")].append(int(leftRepRight[1])/len(currentRef[0]))
            else:
                dict[(currentKey,"exact")].append(int(leftRepRight[1])/len(currentRef[0]))

        #if mone == 3:
        #    break
    mergedFile.close()
    return (dict,locationsList)

(diction, locations) =  readInMergedFile("/home/gabdank/Documents/STR_Attempt/Simulation2/merged.output")
print "chromosome:position\treference\texactRange\texactLeader\tlowerRange\texactValues\tlowerValues"

#print locations

bedFile = open("/home/gabdank/Documents/STR_Attempt/Simulation2/bedfile.bed","w")

for l in locations:
    ref = diction[(l,"reference")]
    refPrint = "{"+str(ref[0])+"}"+str(ref[1])
    print str(l) + "\t" + str(refPrint)
    if getLeadingValue(diction[(l,"exact")]) != None:
        print diction[(l,"exact")]
        '''
        exactDelta = abs(ref[1]-getLeadingValue(diction[(l,"exact")]))
        if exactDelta>3:
            print '*'+str(l[0])+":"+str(l[1])+"\t"+refPrint+"\t"+str(getRange(diction[(l,"exact")]))+"\t"+str(getLeadingValue(diction[(l,"exact")]))+"\t"\
          +str(getRange(diction[(l,"lower")]))+"\t"+str(sorted(diction[(l,"exact")]))+"\t"+str(sorted(diction[(l,"lower")]))
            bedFile.write(str(l[0])[3:]+"\t"+str(l[1])+"\t"+str(int(l[1])+int(ref[1]))+"\n")
        #else:
        #    print str(l[0])+":"+str(l[1])+"\t"+refPrint+"\t"+str(getRange(diction[(l,"exact")]))+"\t"+str(getLeadingValue(diction[(l,"exact")]))+"\t"\
        #      +str(getRange(diction[(l,"lower")]))+"\t"+str(sorted(diction[(l,"exact")]))+"\t"+str(sorted(diction[(l,"lower")]))
        '''
    if getLeadingValue(diction[(l,"lower")]) != None:
        print diction[(l,"lower")]
    '''else:

        lowerUp = -1
        lowerRange = getRange(diction[(l,"lower")])
        if len(lowerRange)==2:
            lowerUp = lowerRange[1]
            if lowerUp-ref[1]>3:
                print '*'+str(l[0])+":"+str(l[1])+"\t"+refPrint+"\t"+str(getRange(diction[(l,"exact")]))+"\t"+str(getLeadingValue(diction[(l,"exact")]))+"\t"\
          +str(getRange(diction[(l,"lower")]))+"\t"+str(sorted(diction[(l,"exact")]))+"\t"+str(sorted(diction[(l,"lower")]))
                bedFile.write(str(l[0])[3:]+"\t"+str(l[1])+"\t"+str(int(l[1])+int(ref[1]))+"\n")
            #else:
            #    print str(l[0])+":"+str(l[1])+"\t"+refPrint+"\t"+str(getRange(diction[(l,"exact")]))+"\t"+str(getLeadingValue(diction[(l,"exact")]))+"\t"\
            #  +str(getRange(diction[(l,"lower")]))+"\t"+str(sorted(diction[(l,"exact")]))+"\t"+str(sorted(diction[(l,"lower")]))
        #else:
        #    print str(l[0])+":"+str(l[1])+"\t"+refPrint+"\t"+str(getRange(diction[(l,"exact")]))+"\t"+str(getLeadingValue(diction[(l,"exact")]))+"\t"\
        #      +str(getRange(diction[(l,"lower")]))+"\t"+str(sorted(diction[(l,"exact")]))+"\t"+str(sorted(diction[(l,"lower")]))
    '''
bedFile.close()