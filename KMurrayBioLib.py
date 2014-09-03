#!/usr/bin/env python
'''
Keith Murray's Bioinformatics Library 
This is intended to become my main import, so this can use and
hold almost every file I need

'''
import math
import numpy
import sys
import os.path


def dnaCharToInt(seq):
    '''
    Converts string of letters (ATCG) into an
    array of numbers and returns it as "dataConverted"
        Base A becomes 0
        Base C becomes 1
        Base G becomes 2
        Base T becomes 3
    '''
    seq = list(seq)
    for i in range(len(seq)):
	if (seq[i] == 'A'):
	    seq[i] = 0
	elif (seq[i] == 'a'):
	    seq[i] = 0
	elif (seq[i] == 'C'):
	    seq[i] = 1
	elif (seq[i] == 'c'):
	    seq[i] = 1
	elif (seq[i] == 'G'):
	    seq[i] = 2
	elif (seq[i] == 'g'):
	    seq[i] = 2
	elif (seq[i] == 'T'):
	    seq[i] = 3
	elif (seq[i] == 't'):
	    seq[i] = 3
	else:
	    seq[i] = 4
    
    return seq


def amiMatrixProbability(intSeq) :
    '''
    Taking in the genetic sequence with bases represented
    as int values, this function outputs an array where each
    value, 0 through 4, represents the probability that a base,
    A, C, G, T, or N will occur
    '''
    fullBaseList = intSeq
    length = len(fullBaseList)
    indexLoopCount = 0
    independentProb = [[0]*5]

    while indexLoopCount < length :
        if fullBaseList[indexLoopCount] == 0 :
            independentProb[0][0] += 1
        elif fullBaseList[indexLoopCount] == 1 :
            independentProb[0][1] += 1
        elif fullBaseList[indexLoopCount] == 2 :
            independentProb[0][2] += 1
        elif fullBaseList[indexLoopCount] == 3 :
            independentProb[0][3] += 1            
        indexLoopCount +=1
    
    i = 0
    a = 0
    for i in range (5) :
        a = independentProb[0][i]
        if ((length != 0)): 
            a = float(a) / length
        independentProb[0][i] = a
        a = 0
        i += 1

    return independentProb
    

def getAMIvector(seq, window) :
    '''
    Summary of the AMI can be found in the following paper
    http://www.biomedcentral.com/1471-2105/9/48
    The Average Mutual Information (AMI) profile is a vector 
    representation of a genetic sequence. Each entry i in the
    vector represents the probability that any base i + 1 
    forward in the sequence can be predicted give that you know
    the current base in the sequence. 
    This function requires the input sequence, and vector size,
    also called the window size. The window size is recommended 
    to be between 64 to 200.
    '''
    intSeq = dnaCharToInt(seq)
    fullBaseList = intSeq
    length = len(fullBaseList)
    independentProb = amiMatrixProbability(intSeq)
    window = window - 1
    k = 0
    amiProfileList = []

    for k in range(1, int(window+2)) :
        i = 0
        j = 0
        jointProb =  [ [ 0. for a in range(5) ] for b in range(5) ]        
        for window in range(length - k ) :
            i = fullBaseList[window]
            j = fullBaseList[window + k]
            jointProb[i][j] += 1
        i = 0
        j = 0
        for i in range(4) :
            for j in range(4) :
                if ((length-k) != 0) :
                    jointProb[i][j] = (jointProb[i][j]) / (length - k) 
        iOfK = 0
        logNum = 0
        for i in range(4) :
            for j in range(4) :
		logNum = 0
                if ((independentProb[0][i] != 0)and(independentProb[0][j] !=0)) : 
                    logNum = (float(jointProb[i][j]) / (float(independentProb[0][i]) * independentProb[0][j]) )  
                if (logNum != 0) :
                    iOfK += (float(jointProb[i][j]) * math.log(logNum))  
        amiProfileList.append(iOfK)   
    return amiProfileList


def vectDistCorrelation(a, b):
    '''

    '''
    i = 0
    j = 0
    meanA = 0
    meanB = 0
    for i in range(len(a)) :
        meanA = meanA + float(a[i])
    for j in range(len(b)) :
        meanB = meanB + float(b[j])
    meanA = float(meanA)/len(a)
    meanB = float(meanB)/len(b)    
    k = 0
    l = 0
    A = 0
    Asqrt = 0
    B = 0
    Bsqrt = 0
    dotProduct = 0
    denom = 0
    distance = 0
    formatDistance = 0
    for k in range(len(a)) :
        A = float(a[k]) #- meanA
        B = float(b[k]) #- meanB
        dotProduct = dotProduct + (A - meanA)*(B - meanB)        
        Asqrt = Asqrt + ((A - meanA) * (A - meanA))
        Bsqrt = Bsqrt + ((B - meanB) * (B - meanB)) 
    Asqrt = math.sqrt(Asqrt)
    Bsqrt = math.sqrt(Bsqrt)
    denom = Asqrt * Bsqrt
    distance = float(dotProduct)/denom
    correlationDistance = 1 - distance 
    return correlationDistance

def vectDistEuclidean(a, b) :
    i = 0
    x = 0
    y = 0
    difference = 0
    euclideanDistance = 0
    sumEuc = 0
    for i in range(len(a)) :
        x = float(a[i])
        y = float(b[i])
        
        difference = (x - y)
        sumEuc += difference*difference
    euclideanDistance = math.sqrt(sumEuc)
    return euclideanDistance

def distanceMatrix(fileName, filePath, metric):
    # Readin Block
        #File Name:
    vectorFile = open(str(filePath) + str(fileName), 'r')
    duplicateCheck = 0
    genomeNameList = []
    hardVector = []

    #listName = []
    listVectors = []
    
    while True :
        hardNameOfFile = vectorFile.readline().strip()
        if (hardNameOfFile == ''):
            vectorFile.close()
            break
        hardLength = vectorFile.readline().strip()
        hardWindow = vectorFile.readline().strip()
        hardVectorRaw = vectorFile.readline().strip()
        hardVectorRaw = hardVectorRaw[1:-1]
        
        hardVector = hardVectorRaw.split(',')
        i = 0
        for i in range(len(hardVector)) :
            hardVector[i].strip()

        genomeNameList.append(hardNameOfFile)

        listVectors.append(hardVector)
        
       
    vectorFile.close()
    distanceMatrix =  [ [ 0. for a in range(len(listVectors)) ] for b in range(len(listVectors)) ]
    #Grab two vectors Block
    i = 0
    j = 0
    distanceXY = 0

    print ("    " + str(len(listVectors)) + " Entries")
    print('Computing Distance Matrix...')

    x = 0


    for i in range(len(listVectors) ) :
        for j in range(len(listVectors)) :#- i -1) :
            #Switch between the two metrics
            if ( (str(metric) == str('c')) or (str(metric) == str('C')) ):
                #print "Correlation Matrix"
                distanceXY = Correlation(listVectors[i], listVectors[j])
            if ( (str(metric) == str('e')) or (str(metric) == str('E')) ):
                #print "Euclidean Matrix"
                distanceXY = euclidean(listVectors[i], listVectors[j])
            distanceMatrix[i][j] = distanceXY
            distanceMatrix[j][i] = distanceXY
        
            
    


    #Save and Output to file block

    i = 0
    j = 0
    if ( (str(metric) == str('c')) or (str(metric) == str('C')) ):
        #print "Correlation Matrix"
	distanceFile = open(str(filePath) + 'AMI_OUTFILE/DistanceMatrixC.txt', "w")
	displayFile = open(str(filePath) + 'AMI_OUTFILE/DisplayFileC.txt', "w")
	labelFile = open(str(filePath) + 'AMI_OUTFILE/LabelFileC.txt', "w")
    if ( (str(metric) == str('e')) or (str(metric) == str('E')) ): 
	distanceFile = open(str(filePath) + 'AMI_OUTFILE/DistanceMatrixE.txt', "w")
	displayFile = open(str(filePath) + 'AMI_OUTFILE/DisplayFileE.txt', "w")
	labelFile = open(str(filePath) + 'AMI_OUTFILE/LabelFileE.txt', "w")
    
    distanceFile.write(str(len(listVectors)))
    distanceFile.write("\n")

    loopTrue = 0
    nameArray = []
    
    
    for i in range(len(listVectors) ) :
        name = 'seq' + str(i)
        distanceFile.write(str(name))
        distanceFile.write("    ")
        for j in range(len(listVectors) ) :
            distanceFile.write(str(distanceMatrix[i][j]))
            distanceFile.write("    ")
        distanceFile.write("\n")
        if (genomeNameList[i][0] == 'g'):
            title = 1
        elif (genomeNameList[i][0] == '_'):
            title = 2
	else :
	    #THIS ELSE STATEMENT WAS A MODIFICATION FOR DICTIONARIES: 
	    #               FIX AND REMOVE LATER
	    title = 0
       
        match = 0
        tempHold = []
        genomeName = genomeNameList[i]
        tempHold = genomeName.split(' ')
	#print title
	#print tempHold
        location = 0
        k = 0
        if (loopTrue == 0) :
            nameArray.append(tempHold[title])

        for k in range (len(nameArray)) :
            if (tempHold[title] == nameArray[k]) :
                match = 1
                location = k

        if (match == 1) :
            line = str(location) + " 0 4\n"
        elif (match == 0) :
            nameArray.append(tempHold[title])
            line = str((len(nameArray) - 1)) + " 0 4\n"
        else :
            print "CRAKE WAS HERE and you have an error"
            

        loopTrue +=1
        displayFile.write(line)
        #print(genomeNameList[i][0])
        labelFile.write(genomeNameList[i])
        labelFile.write("\n")

    distanceFile.close()
    print("Matrix Done.")
    return



























def help(function =""):
    
    if ((function == "dnaCharToInt") or (function == 0)):
	msg = '''
    *******************************************
    FUNCTION HELP:	0:  dnaCharToInt
    IN:			DNA sequence (String/List)
    OUT: 		ints represent DNA string (List)

    TIME COMPLEXITY: 	O(n)
			n = sequence length

    SUMMARY:
    Converts string of letters (ATCG) into an
    array of numbers and returns it as "dataConverted"
        Base A becomes 0
        Base T becomes 1
        Base G becomes 2
        Base C becomes 3

    *******************************************
	
	'''
    elif ((function == "amiMatrixProbability") or (function == 1)):
	msg = '''
    *******************************************
    FUNCTION HELP:	1:  amiMatrixProbability
    IN:			Sequence as array of ints [List of Ints]
    OUT:		Independent Probability [List ]

    TIME COMPLEXITY:	O(n)
			n = sequence length

    SUMMARY:
    Taking in the genetic sequence with bases represented
    as int values, this function outputs an array where each
    value, 0 through 4, represents the probability that a base,
    A, C, G, T, or N will occur.

    *******************************************
	    '''

    elif ((function == "getAMIvector") or (function == 2)):
	msg = '''
    *******************************************
    FUNCTION HELP:	2:  getAMIvector
    IN:			Sequence [String]
			Window length [Int]
    OUT:		AMI Vector [List of Floats]

    TIME COMPLEXITY:	O(nw)
			n = sequence length
			w = window length

    SUMMARY:
    Summary of the AMI can be found in the following paper
    http://www.biomedcentral.com/1471-2105/9/48
    The Average Mutual Information (AMI) profile is a vector 
    representation of a genetic sequence. Each entry i in the
    vector represents the probability that any base i + 1 
    forward in the sequence can be predicted give that you know
    the current base in the sequence. 
    This function requires the input sequence, and vector size,
    also called the window size. The window size is recommended 
    to be between 64 to 200.
    
    *******************************************
	    '''

	
    else:
	msg = '''
    *******************************************
      ***Please enter one of the following***
      ***  functions for more information *** 
0:  dnaCharToInt
1:  amiMatrixProbability
2:  getAMIvector


	'''

    print msg
    return


