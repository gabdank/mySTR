'''

First read in all the information.
You will need for the further analysis for any given STR genomic location, be able to produce:
1. chromosome:posiotion
2. For any given strain (un-dependent on the "interest" the list of repetitive region length and the alignments.
3. Ability to summarize the alignments repetitive region length.


'''
from operator import itemgetter


def initiateReferenceDictionary(intDict, listOfFileNames, flankLengthThreshold):
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

def processRecord(record, flanksLengthThreshold):
    ar = record.split("${")
    referenceData = ar[0]
    arra = referenceData.split("$")
    chromo = arra[5].split(":")[1]
    position = arra[6].split(":")[1]
    repetitiveSequenceLength = getReferenceRepeatLength(referenceData)
    unitLength = getReferenceUnitLength(referenceData)
    fullRepetitiveSequenceLength = (int(repetitiveSequenceLength)/int(unitLength))*unitLength
    delta = repetitiveSequenceLength-fullRepetitiveSequenceLength
    return (chromo,position,processReferenceData(referenceData, flanksLengthThreshold, delta))


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

def processNumbers(strainFile, locationsDictionary, strainName, strainSpecificDiction):
    m = 0
    currentRecord = ""
    fil = open(strainFile, "r")
    for l in fil:
        m+=1
        n = 0
        if l[0:14]=="repetitiveUnit":
            if len(currentRecord)>0:
                locationsDictionary[getKey(currentRecord[0:-1])][strainName]=getNumbers(currentRecord[0:-1])
                refNumOfReps = locationsDictionary[getKey(currentRecord[0:-1])]['reference'][2]
                listExacts = locationsDictionary[getKey(currentRecord[0:-1])][strainName][0]
                listLowers = locationsDictionary[getKey(currentRecord[0:-1])][strainName][1]

                if isInteresting(listExacts, listLowers, refNumOfReps)== True:
                    strainSpecificDiction[getKey(currentRecord[0:-1])]=currentRecord

                currentRecord = l.strip()
            else:
                currentRecord = l.strip()+"$"
        else:
            currentRecord+= l.strip()+"$"
    locationsDictionary[getKey(currentRecord)][strainName]=getNumbers(currentRecord)
    refNumOfReps = locationsDictionary[getKey(currentRecord)]['reference'][2]
    listExacts =  locationsDictionary[getKey(currentRecord)][strainName][0]
    listLowers =  locationsDictionary[getKey(currentRecord)][strainName][1]
    if isInteresting(listExacts, listLowers, refNumOfReps)== True:
        strainSpecificDiction[getKey(currentRecord)]=currentRecord

    fil.close()

def isInteresting(exactValues, lowerBoundValues, referenceValue):
    # check if there are at least 3 values different from the reference in the exact values
    # find the majority by doing a small histo
    valDic = {}
    for v in exactValues:
        if int(v) in valDic:
            valDic[int(v)] += 1
        else:
            valDic[int(v)] = 1
    max = 0
    for v in valDic:
        if valDic[v]>max:
            max = valDic[v]

    maxCounter = 0
    location = 0
    for v in valDic:
        if max - valDic[v] < 2:
            maxCounter += 1
            location = v
    if maxCounter == 1 and valDic[location]>3:
        if abs(location-int(referenceValue))>3:
            return True
        if abs(location-int(referenceValue))<2 and valDic[location]>5: # updated version that ensures we will not have multiple (>5) reference supportive reads, but only different from it
            return False
    mone = 0
    for v in lowerBoundValues:
        if int(v)>(int(referenceValue)+2):
            mone += 1
    if mone >= 3:
        return True
    return False

def extractCommonLength(exactValues, lowerBoundValues, referenceValue): #like isInteresting - will attempt to extract th eprobable value of the STR expansion/contraction, from the exact and lower bound values
    valDic = {}
    for v in exactValues:
        if int(v) in valDic:
            valDic[int(v)] += 1
        else:
            valDic[int(v)] = 1
    max = 0
    for v in valDic:
        if valDic[v]>max:
            max = valDic[v]

    maxCounter = 0
    location = 0
    for v in valDic:
        if max - valDic[v] < 2:
            maxCounter += 1
            location = v
    return location








genomicLocationsDictionary = {}

filesList = ["/media/gabdank/Disk3/MiliionMutation/N2/AlignmentData.results", "/media/gabdank/Disk3/MiliionMutation/MY1/AlignmentData.results", "/media/gabdank/Disk3/MiliionMutation/MY2/AlignmentData.results"]

# Initiating the genomic locations dictionary, that will contain refrence genome information for every location
# the dictionary will contain all locations that were found to be interesting in different strain runs. So it is ESSENTIAL to run all the strains results to initiate the reference list,
# otherwise there would be missing locations that will cause exceptions later.
# each entry in the dictionary looks like that: ('II', '12309790'): {'reference': ('TTTTTTTATACGATTTTTTA', 'AAT', 9, 'ZERO', 'TTATTACGTTTGCATTGAAA')}

initiateReferenceDictionary(genomicLocationsDictionary, filesList, 20)


# initiation of special dictionaries - that will store the "interesting locations data"
# when they would be constructed the genomicLocationsReferenceDictionary will be updated as well.
# It will sotre the exact and lower bound repeat lengths (if there were any reads alignable to this place) - that will allow further
# processing. The special dictionaries could be used for the "alignments" - the numbers, stored in
# genomicLocationsReferenceDictionary could be used for preparation of the TABLE of the results.

strainDictions = {}
strainDictions['N2']={}
strainDictions['MY1']={}
strainDictions['MY2']={}

strainFileName = "/media/gabdank/Disk3/MiliionMutation/N2/AlignmentData.results"
processNumbers(strainFileName, genomicLocationsDictionary,"N2",strainDictions['N2'])

strainFileName = "/media/gabdank/Disk3/MiliionMutation/MY1/AlignmentData.results"
processNumbers(strainFileName, genomicLocationsDictionary,"MY1",strainDictions['MY1'])

strainFileName = "/media/gabdank/Disk3/MiliionMutation/MY2/AlignmentData.results"
processNumbers(strainFileName, genomicLocationsDictionary,"MY2",strainDictions['MY2'])

# You have to remember that dictionary may or may not contain info for given strain for given location
# The only thing promised is the presence of the reference genome STR description.
# Other than that you have to make sure there is data (in terms of alignments) at all






#####$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

def printAlignmentData(data, unitLength, threshold, referenceLength, referenceRepeat):
    #print "PRINTING ALIGNMENT"

    #getRepeatInfo(data, unitLength, threshold, referenceLength, referenceRepeat)

    #print "END OF THE REPAIR"


    array = data.split("$")
    if len(array)>2:
        left = array[0].split(":")[1]
        right = array[1].split(":")[1]
        repeat = array[2].split(":")[1]
        repeatUnit = array[2].split(":")[1][0:unitLength]
        readNumOfRepeats = len(repeat)/unitLength
        referenceNumberOfRepeats = referenceLength/unitLength

        absoluteDifference = abs(referenceNumberOfRepeats-readNumOfRepeats)

        if repeatUnit==referenceRepeat:
            lF = leftFlankTreatment(left, threshold)
            rF = rightFlankTreatment(right, threshold)

            if (len(left)<10 or len(right)<10) and readNumOfRepeats>referenceNumberOfRepeats+2 :
                print lF + "\t{"+repeatUnit+"}"+str(readNumOfRepeats)+"\t"+rF
            else:
                if len(left)>=10 and len(right)>=10 and absoluteDifference>1:

                    print lF + "\t{"+repeatUnit+"}"+str(readNumOfRepeats)+"\t"+rF
        else:
            tempRepeat = repeatUnit
            tzuza = 0
            flag = False
            repeatProblem = False

            while flag!=True:
                tzuza += 1
                tempRepeat = tempRepeat[1:]+tempRepeat[:1]
                if tempRepeat==referenceRepeat:
                    flag=True
                if tzuza==len(tempRepeat):
                    flag = True
                    repeatProblem = True
            if repeatProblem == False:
                ostatok = len(repeatUnit)-tzuza
                possibleRepeat = left[-ostatok:]+repeat[:tzuza]
                updatedLeftFlank = ""
                combinedRepeat = "zopa"

                if possibleRepeat == referenceRepeat:
                    combinedRepeat = possibleRepeat+repeat[tzuza:]
                    updatedLeftFlank = left[:-ostatok]

                else:
                    combinedRepeat = repeat[tzuza:]
                    updatedLeftFlank = left+repeat[:tzuza]

                readNumOfRepeats = len(combinedRepeat)/len(repeatUnit)

                lF = leftFlankTreatment(updatedLeftFlank, threshold)
                rF = rightFlankTreatment(right, threshold)
                absoluteDifference = abs(referenceNumberOfRepeats-readNumOfRepeats)
                if (len(left)<10 or len(right)<10) and readNumOfRepeats>referenceNumberOfRepeats+2:
                    print lF + "\t{"+referenceRepeat+"}"+str(readNumOfRepeats)+"\t"+rF
                else:
                    if len(left)>=10 and len(right)>=10 and absoluteDifference>1:
                        print lF + "\t{"+referenceRepeat+"}"+str(readNumOfRepeats)+"\t"+rF



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



#########################$$$$$$$$$$

def printBeauty(t):
    return t[0]+"\t{"+t[1]+"}"+str(t[2])+"\t{"+t[1]+"}"+str(t[3])+"\t{"+t[1]+"}"+str(t[4])+"\t{"+t[1]+"}"+str(t[5])




moneOfPrinted= 0
n2_counter = 0
my1_counter = 0

listOfStrains = ['N2','MY1','MY2']

listOfTuples = []

for k in genomicLocationsDictionary:
    d = genomicLocationsDictionary[k]

    strains = d.keys()
    refNumOfReps = d['reference'][2]
    #print d

    #print "Reference number of reps:"+str(refNumOfReps)
    lengthReport = {}
    for strainName in listOfStrains:
        lengthReport[strainName]="-"
    for strainName in listOfStrains:
        if strainName in d:
            if len(d[strainName][0])>2:
                exacts = d[strainName][0]
                lowers = d[strainName][1]
                lengthReport[strainName] = str(extractCommonLength(exacts,lowers,refNumOfReps))

    interestingCounter = 0
    printFlag = False
    for strain in strains: # walking over the strains for this location

        if strain != "reference":
            exacts = d[strain][0]
            lowers = d[strain][1]


            if isInteresting(exacts,lowers,refNumOfReps)==True:
                locDa = k
                reDa = genomicLocationsDictionary[k]['reference']
                refPrint = "chr"+str(locDa[0])+":"+str(locDa[1])+"\t{"+str(reDa[1])+"}"+str(reDa[2])
                refPrint += "\t{"+str(reDa[1])+"}"+str(lengthReport['N2']) +"\t{"+str(reDa[1])+"}"+str(lengthReport['MY1'])
                refPrint += "\t{"+str(reDa[1])+"}"+str(lengthReport['MY2'])
                locus = "chr"+str(locDa[0])+":"+str(locDa[1])
                rUnit = reDa[1]
                if lengthReport['N2']!= '-':
                    n2Delta = abs(int(lengthReport['N2'])-int(reDa[2]))
                else:
                    n2Delta = 0

                if lengthReport['MY1']!= '-':
                    my1Delta = abs(int(lengthReport['MY1'])-int(reDa[2]))
                else:
                    my1Delta = 0

                if lengthReport['MY2']!= '-':
                    my2Delta = abs(int(lengthReport['MY2'])-int(reDa[2]))
                else:
                    my2Delta = 0

                tupla = (locus,rUnit,reDa[2],lengthReport['N2'],lengthReport['MY1'], lengthReport['MY2'], n2Delta, my1Delta, my2Delta , max( n2Delta, my1Delta, my2Delta),d)

                listOfTuples.append(tupla)

                interestingCounter += 1


    if True: #interestingCounter>0 and printFlag == True:
        moneOfPrinted += 1
        #print "~~~~~~~~~~~~~~~~~~"
        #print refPrint

        for strain in strains:
            if strain != "reference":
                exacts = d[strain][0]
                lowers = d[strain][1]
                if isInteresting(exacts,lowers,refNumOfReps)==True:


                    locationData = k
                    refData = genomicLocationsDictionary[k]['reference']

                    print "EXTRACTED LENGTH:\t{"+refData[1]+"}"+str(extractCommonLength(exacts,lowers,refNumOfReps))
                    print "\nStrain:"+strain+"\tExact:"+str(genomicLocationsDictionary[k][strain][0])+ " Lower:"+str(genomicLocationsDictionary[k][strain][1])+"\n"

                    if (refData[3]=='ZERO'):
                        print refData[0]+"\t{"+refData[1]+"}"+str(refData[2])+"\t"+refData[4]
                    else:
                        print refData[0]+"\t{"+refData[1]+"}"+str(refData[2])+"\t"+refData[4]+"\t"+"Tail="+refData[3]
                    print "--------------------------------------------------------------"

                    ar =  strainDictions[strain][k].split("${")
                    referenceData = ar[0]
                    alignmentsData = ar[1][:-1]

                    for element in  alignmentsData.split("]#"):
                        if (len(element.strip())>0):
                            printAlignmentData(element[1:],getReferenceUnitLength(referenceData), 20, getReferenceRepeatLength(referenceData),genomicLocationsDictionary[k]['reference'][1] )
        print "==========================================\n"

print "chromosome:position\tReference\tN2\tMY1\tMY2"
x = sorted(listOfTuples, key=itemgetter(9))
for ele in reversed(x):
    #print ele
    print printBeauty(ele)
    #print ele[10]