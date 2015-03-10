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

   
def processNumbers(strainFile, locationsDictionary, strainName, strainSpecificDiction):
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
            print currentRecordKey
            #print locationsDictionary[currentRecordKey]
        #if l[0]=='{': #starting section of aligned reads
        #if l[0]=='}': #finished to read section of aligned reads
        if l[0]=="[":
            print l[1:-2].split('$')
            
        '''m+=1
        n = 0
        if l[0:14]=="repetitiveUnit":
            if len(currentRecord)>0:
                locationsDictionary[getKey(currentRecord[0:-1])][strainName]=getNumbers(currentRecord[0:-1])
                refNumOfReps = locationsDictionary[getKey(currentRecord[0:-1])]['reference'][2]
                listExacts = locationsDictionary[getKey(currentRecord[0:-1])][strainName][0]
                listLowers = locationsDictionary[getKey(currentRecord[0:-1])][strainName][1]

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

    
dataDictionary = {}
listOfMergedFiles = ['merged.output','merged.output']

initRefDict(dataDictionary,listOfMergedFiles,20)

# at this stage we have finished to extract all the reference information and now we should extract strain specific information

for entry in dataDictionary:
    print entry
    print dataDictionary[entry]

strainDictionary = {}
strainDictionary["sampleStrain"]={}
strainFileName = "merged.output"
processNumbers(strainFileName, dataDictionary, "sampleStrain", strainDictionary["sampleStrain"])

    
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