# logical plan of operation
# 1. Read the merged file and create a dictionary with all reference repeat lengths and information
# After parsing the merged file you should be able to get from the dictionary everything usingn chromosome 
# and position as a key.
# For reference what would be stored is: 
# (chr,position, (left flank, repeatUnit, length, 'ZERO' or number of partial unit , right flank))
# given number of libraries to process we have to create a dictionary of reference genome for all of 
# them, even if some doesn't contain information of the position - due to the lack of reads.


''' function goes over ALL detection files (merged type files) and  builds a reference dictionary
allowing down the road to extract information about different reference locations
it stores for every location (chr,position) a dictioanry with key that correspond to the strains
names. For each location we can retreive information of the refrence genome or any given strain, but
the function will only build the refrence currently.
'''

import numpy as np

def initRefDict(diction, fileNamesList, flankLengthThreshold):    
    for name in fileNamesList:
        mergedFile = open(name,"r")        
        print "Extracting from "+name
        for line in mergedFile:           
            if line[0:1]=='r':               
                extrData = extractReferenceInfoFromRecord(line, flankLengthThreshold)           
                key = (extrData[0],extrData[1])
                dictionInternal = {}
                dictionInternal['reference'] = extrData[2]
                if not key in diction:
                    diction[key]=dictionInternal                       
        mergedFile.close()
    #print diction

def extractReferenceInfoFromRecord(record, flanksLengthThreshold):
    arra = record.split("$")
    chromo = arra[5].split(":")[1]
    position = arra[6].split(":")[1]
    repetitiveSequenceLength = getReferenceRepeatLength(record)
    unitLength = getReferenceUnitLength(record)
    fullRepetitiveSequenceLength = (int(repetitiveSequenceLength)/int(unitLength))*unitLength
    delta = repetitiveSequenceLength-fullRepetitiveSequenceLength
    return (chromo,position,processReferenceData(record, flanksLengthThreshold, delta))
   

def getReferenceUnitLength(referenceData):
    array = referenceData.split("$")
    return len(array[0].split(":")[1])

def getReferenceRepeatLength(referenceData):
    array = referenceData.split("$")
    return len(array[4].split(":")[1])

def getRepetitiveSequence(referenceData):
    array = referenceData.split("$")
    return array[4].split(":")[1]

def processReferenceData(data, threshold, delta):
    array = data.split("$")
    unit = array[0].split(":")[1]
    repLength = int(array[7].split(":")[1])
    numOfReps = repLength/len(unit)
    actualRepeat =  array[4].split(":")[1][0:len(unit)]

    unfinishedBusiness = "ZERO"
    if delta>0:
        unfinishedBusiness = getRepetitiveSequence(data)[-delta:]

    return (array[2].split(":")[1][-threshold:], actualRepeat, numOfReps ,unfinishedBusiness, array[3].split(":")[1][0:threshold])

   
def processNumbers(strainFile, locationsDictionary, strainName):
    m = 0
    currentRecordKey = ""
    fil = open(strainFile, "r")
    #The expectation is to get file of the following format:
    '''
    repetitiveUnit:
    {
    [leftSequence:XXX$...]
    []
    []
    []
    }
    repetitiveUnit:    
    
    
    '''
    for l in fil:
        if l[0]=='r': # we are dealing with beginning of location
            array = l.strip()[:-1].split('$')
            #print array
            #print array[-3]
            #print array[-2]
            currentRecordKey = (array[-3].split(":")[1],array[-2].split(":")[1])

            
            #print ""
            #print currentRecordKey
            #print locationsDictionary[currentRecordKey]
            #print ""
            genomicLeft = locationsDictionary[currentRecordKey]['reference'][0]
            genomicRight = locationsDictionary[currentRecordKey]['reference'][4]
            genomicUnit = locationsDictionary[currentRecordKey]['reference'][1]
            genomicSheerit = locationsDictionary[currentRecordKey]['reference'][3]
            genomicNumberOfReps = int(locationsDictionary[currentRecordKey]['reference'][2])
            genomicRepetitive = ''
            mona = 0            
            while mona<genomicNumberOfReps:
                mona += 1
                genomicRepetitive += genomicUnit
            if genomicSheerit!='ZERO':
                genomicRepetitive +=genomicSheerit
                
                
            #print "genomic repetative = "+genomicRepetitive
            #print "genomic construction:"
            #print genomicLeft + "\t" +genomicRepetitive[0:8]+"..."+genomicRepetitive[-8:]+"\t"+str(len(genomicRepetitive))+"\t"+genomicRight
            #print "------------------------------------------------------------------------------------------"
            # initiate two lists - for exact and for lower values:
            exactList = []
            lowerList = []            
        if l[0]=="[":
            #print l[1:-2].split('$')
            localArr = l[1:-2].split('$')
            localLeft = localArr[0].split(':')[1]
            localRight = localArr[1].split(':')[1]
            localRepetitive = localArr[2].split(':')[1]
            
            #print "local repetitive : "+localRepetitive
            #print localLeft+"\t"+localRepetitive[0:8]+"..."+localRepetitive[-8:]+"\t"+str(len(localRepetitive))+"\t"+localRight
            localUnit = localRepetitive[0:len(genomicUnit)]
            
            leftIndex = leftAlignment(genomicLeft, localLeft)
            #print leftIndex
            rightIndex = rightAlignment(genomicRight, localRight)
            #print rightIndex
            if leftIndex ==0:
                if len(localLeft)<=20:
                    storedLeft = localLeft
                else:
                    storedLeft = localLeft[-20:]
                combinedLocalLeft = leftFlankTreatment(localLeft,20)
                combinedLocalRepetitive = localRepetitive+localRight[:rightIndex] 
            else:          
                shortenedLeft = localLeft[:-leftIndex]
                if len(shortenedLeft)<=20:
                    storedLeft = shortenedLeft
                else:
                    storedLeft = shortenedLeft[-20:]
                combinedLocalLeft = leftFlankTreatment(shortenedLeft,20) 
                combinedLocalRepetitive = localLeft[-leftIndex:]+localRepetitive+localRight[:rightIndex] 
            shortenedRight = localRight[rightIndex:]
            if len(shortenedRight)<=20:
                storedRight = shortenedRight
            else:
                storedRight = shortenedRight[:20]
            recordToStore = (storedLeft, combinedLocalRepetitive, storedRight)
            
            combinedLocalRight = rightFlankTreatment(localRight[rightIndex:],20)       
            combinedRepeatLength = len(combinedLocalRepetitive)
            #print combinedLocalLeft+"\t"+combinedLocalRepetitive[:8]+"..."+combinedLocalRepetitive[-8:]+"\t"+str(len(combinedLocalRepetitive))+"\t"+combinedLocalRight
            if len(storedLeft)<8 or len(storedRight)<8: #lower case
                lowerList.append(recordToStore)                
            else:
                exactList.append(recordToStore)
            #print localAlignment(genomicLeft,genomicRepetitive,genomicRight, localLeft,localRepetitive, localRight)            
            #print leftIndex
            #print genomicLeft
            #print leftFlankTreatment(localLeft[:-leftIndex],20)            
            #print genomicRepetitive
            #print combinedLocalRepetitive
            #print rightIndex
            #print genomicRight
            #print rightFlankTreatment(localRight[rightIndex:],20)
            #break
            #print localUnit
            #print localRepetitive            
            #print localRight
            #print localLeft
            #print "=============================================="
            
            
        if l[0]=="}": #time to put all the data in the aligned sequences into dictioanry entry
            locationsDictionary[currentRecordKey][strainName]=(exactList,lowerList)
            
            '''exactLengthsList = extractLengthsFromValuesList(exactList)
            lowerLengthsList = extractLengthsFromValuesList(lowerList)
            exactArray = np.array(exactLengthsList)            
            lowerArray = np.array(lowerLengthsList)            
            
            print exactArray
            print np.average(exactArray)
            print np.std(exactArray)
            print lowerArray
            print np.average(lowerArray)
            print np.std(lowerArray)
            '''
        #if l[0]=='{': #starting section of aligned reads
        #if l[0]=='}': #finished to read section of aligned reads
        '''if l[0]=="[":
            print l[1:-2].split('$')
            
        m+=1
        n = 0
        if l[0:14]=="repetitiveUnit":
            if len(currentRecord)>0:
                locationsDictionary[getKey(currentRecord[0:-1])][strainName]=getNumbers(currentRecord[0:-1])
                refNumOfReps = locationsDictionary[getKey(currentRecord[0:-1])]['reference'][2]
                listExacts = locationsDictionary[getKey(currentRecord[0:-1])][strainName][0]
                listLowers = locationsDictionary[getKey(currentRecord[0:-1])][strainName][1]
                #if isInteresting(listExacts, listLowers, refNumOfReps)== True:
                #    strainSpecificDiction[getKey(currentRecord[0:-1])]=currentRecord


                #if isInteresting(listExacts, listLowers, refNumOfReps)== True:
                #    strainSpecificDiction[getKey(currentRecord[0:-1])]=currentRecord

                currentRecord = l.strip()
            else:
                currentRecord = l.strip()+"$"
        else:
            currentRecord+= l.strip()+"$"
        '''
    #locationsDictionary[getKey(currentRecord)][strainName]=getNumbers(currentRecord)
    #refNumOfReps = locationsDictionary[getKey(currentRecord)]['reference'][2]
    #listExacts =  locationsDictionary[getKey(currentRecord)][strainName][0]
    #listLowers =  locationsDictionary[getKey(currentRecord)][strainName][1]
   
    #if isInteresting(listExacts, listLowers, refNumOfReps)== True:
    #    strainSpecificDiction[getKey(currentRecord)]=currentRecord

    fil.close()


def extractLengthsFromValuesList(anyList):
    vals = []
    for (a,b,c) in anyList:
        vals.append(len(b))
    #print vals
    return vals
    

def leftAlignment(gLeft, lLeft):
    if len(lLeft)>20:
        lLeft = lLeft[-20:]
    scores = {}
    leftShortening = 0
    while leftShortening<10:
        x = -1
        mone = 0
        while x>=-len(lLeft):
            if lLeft[x]==gLeft[x]:
                mone += 1
            x -= 1
        scores[leftShortening]=mone
        lLeft = lLeft[:-1]
        leftShortening += 1

    maxIndex = 0
    maxScore = scores[0]
    for a in scores:
        if scores[a]>maxScore:
            maxScore = scores[a]
            maxIndex = a
    return maxIndex

def rightAlignment(gRight, lRight):
    if len(lRight)>20:
        lRight = lRight[0:20]
    scores = {}
    rightShortening = 0
    while rightShortening<10:
        x = 0
        mone = 0
        while x < len(lRight):
            ##print lRight
            #print x
            if lRight[x]==gRight[x]:
                mone += 1
            x += 1
        scores[rightShortening]=mone
        lRight = lRight[1:]
        rightShortening += 1
    #print scores
    maxIndex = 0
    maxScore = scores[0]
    for a in scores:
        if scores[a]>maxScore:
            maxScore = scores[a]
            maxIndex = a
    return maxIndex


    
def rotateCheck(one, two):
    #print "Word :"+one
    dictOne = {}
    for a in one:
        if not a in dictOne:
            dictOne[a]=1
        else:
            dictOne[a]+=1
    dictTwo = {}
    for b in two:
        if not b in dictTwo:
            dictTwo[b]=1
        else:
            dictTwo[b]+=1
    
    for a in dictOne.keys():
        if a in dictTwo:
            if dictTwo[a]!=dictOne[a]:
                return False
        else:
            return False
    for b in dictTwo.keys():
        if b in dictOne:
            if dictOne[b]!=dictTwo[b]:
                return False
        else:
            return False
            
    #print "-----------"
    return True
    
    
def getKey(record):
    ar = record.split("${")
    referenceData = ar[0]
    arra = referenceData.split("$")
    chromo = arra[5].split(":")[1]
    position = arra[6].split(":")[1]
    return (chromo,position)

def getNumbers(record):
    ar = record.split("${")
    referenceData = ar[0]
    alignmentsData = ar[1]
    arra = referenceData.split("$")
    chromo = arra[5].split(":")[1]
    position = arra[6].split(":")[1]
    repetitiveSequenceLength = getReferenceRepeatLength(referenceData)
    unitLength = getReferenceUnitLength(referenceData)
    fullRepetitiveSequenceLength = (int(repetitiveSequenceLength)/int(unitLength))*unitLength
    delta = repetitiveSequenceLength-fullRepetitiveSequenceLength
    return printNumbers(referenceData, alignmentsData, getReferenceUnitLength(referenceData))

def printNumbers(refD, aliD, unitLength):
    array = refD.split("$")
    unit = array[0].split(":")[1]
    exactVals = []
    lowerVals = []
    for element in  aliD.split("]#"):
        if len(element)>2:
            array = element[1:].split("$")
            left = array[0].split(":")[1]
            right = array[1].split(":")[1]
            numOfRepeats = len(array[2].split(":")[1])/unitLength
            if len(left)<10 or len(right)<10:
                lowerVals.append(numOfRepeats)
            else:
                exactVals.append(numOfRepeats)
    return (sorted(exactVals),sorted(lowerVals))

def leftFlankTreatment(seq, threshold):
    if len(seq)>=threshold:
        return seq[-threshold:]
    else:
        delta = threshold-len(seq)
        toReturn = ""
        for i in range(delta):
            toReturn += "-"
        return toReturn+seq

def rightFlankTreatment(seq,threshold):
    if len(seq)>=threshold:
        return seq[0:threshold]
    else:
        delta = threshold-len(seq)
        toReturn = seq
        for i in range(delta):
            toReturn += "-"
        return toReturn
        
dataDictionary = {}
listOfMergedFiles = ['/home/gabdank/Documents/STR_Attempt/Simulation3/merged.output']

 
initRefDict(dataDictionary,listOfMergedFiles,20)

# at this stage we have finished to extract all the reference information and now we should extract strain specific information
#for entry in dataDictionary:
#    print entry
#    print dataDictionary[entry]

strainFileName = "/home/gabdank/Documents/STR_Attempt/Simulation3/merged.output"
processNumbers(strainFileName, dataDictionary, "sampleStrain")

for (chromo,location) in dataDictionary:
    k = (chromo,location)
    referenceData = dataDictionary[k]['reference']
    if referenceData[3]=='ZERO':
        referenceLength = referenceData[2]*len(referenceData[1])
    else:
        referenceLength = (referenceData[2]*len(referenceData[1]))+len(referenceData[3])

    strainData = dataDictionary[k]['sampleStrain']
    exactL = strainData[0]
    lowerL = strainData[1]

    exactLengthsList = extractLengthsFromValuesList(exactL)
    lowerLengthsList = extractLengthsFromValuesList(lowerL)    
    
    exactArray = np.array(exactLengthsList)            
    #lowerArray = np.array(lowerLengthsList)            
    exactAVG = np.average(exactArray)
    exactSTD = np.std(exactArray)
    #lowerAVG = np.average(lowerArray)
    #lowerSTD = np.std(lowerArray)
    
    #print referenceData
    lineToPrint = chromo +":"+location+"\t["+referenceData[1]+"]"+str(referenceLength)+"\t"
    
    if len(exactLengthsList)>3:
        lineToPrint += str(exactAVG)+":"+str(exactSTD)+"\t"
    else:
        lineToPrint += "-:-\t"
        
    interestingLower = False
    for x in lowerLengthsList:
        if x>referenceLength+8:
            interestingLower=True
    if interestingLower==True:
        lineToPrint += "Long-Lower\t"
    else:
        lineToPrint += "-\t"
    print lineToPrint
    #print exactLengthsList
    #print lowerLengthsList
    #print dataDictionary[(chromo,location)].keys()

    
'''def initiateReferenceDictionary(intDict, listOfFileNames, flankLengthThreshold):
    for name in listOfFileNames:
        alignmentsFile = open(name, "r")
        m = 0
        n = 0
        currentRecord = ""
        for l in alignmentsFile:
            m+=1
            if l[0:14]=="repetitiveUnit":
                if len(currentRecord)>0:
                    answer = processRecord(currentRecord[0:-1], flankLengthThreshold)
                    
                    key = (answer[0],answer[1])
                    print key
                    #print currentRecord
                    # print answer
                    #break
                    diction = {}
                    refData = answer[2]
                    diction["reference"] = refData
                    intDict[key]=diction
                    currentRecord = l.strip()
                else:
                    currentRecord = l.strip()+"$"
            else:
                currentRecord+= l.strip()+"$"


        answer = processRecord(currentRecord, flankLengthThreshold)
        key = (answer[0],answer[1])
        diction = {}
        refData = answer[2]
        diction["reference"] = refData
        intDict[key]=diction
        alignmentsFile.close()
        print intDict

'''
