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

   

dataDictionary = {}
listOfMergedFiles = ['/home/gabdank/Documents/STR_Attempt/Simulation2/merged.output','/home/gabdank/Documents/STR_Attempt/Simulation2/merged.output']
initRefDict(dataDictionary,listOfMergedFiles,20)

for entry in dataDictionary:
    print entry
    print dataDictionary[entry]
    
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