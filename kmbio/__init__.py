#!/usr/bin/env python
import math
from math import sqrt
import os.path
from Bio import SeqIO
import numpy
import scipy
import sys, os
import gzip,cPickle
import re
from subprocess import Popen, PIPE, STDOUT

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import _DistanceMatrix
from ete3 import Tree

'''
kmbio
My cheap personal library for research
======================================

Author(s): Keith Murray

Certain Functions/snipits from:	 Ufuk Nalbantoglu
				 Sam Way

Contact: kmurrayis@gmail.com


*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*
* "Hofstadter's Law: It always takes longer than you expect,           *
*   even when you take into account Hofstadter's Law."                 *
*   --Douglas Hofstadter                                               *
*                                                                      *
* http://xkcd.com/409/                                                 *
*                                                                      *
* I'm really avoiding work right now aren't I?                         *
* ...But I do need to get my charger, so I'll avoid it a bit longer..  *
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*

All of the following code is written for DNA, not amino acids

TODO:
Maybe redo setSequences structure
--Goal:
	-Make developing a new profile easier
Make it simple to see the relationship between a set of sequences 
   ~Around 10 lines of code ideally 


Move code to 3.6 eventually 
--pay specific attention to the use of ordered dictionaries
  --minimize redundancy with order, reoptimize
    ^ This should cut memory usage and improve speed (but speed should be 
	barely noticable) 
--add a 'whodis' function to get a basic readout of a seq or set of seqs
  



'''


class sequence(object):
    # Should really be built on Bio.SeqIO, but that's for later
    def __init__(self, header, seq):
	# sequence needs the seq header, then the seq itself
	self.header = header
	self.seq = seq
	self._setAMI = False
	self._setDBP1 = False
	self._setDBP2 = False
	self._setDict = False
	self._setKmer = False
	self._kVal = 0
    def __getitem__(self, a):
	return self.seq[a]
    def __len__(self):
	return len(self.seq)
    def _dbpload(self, filename):
	file = gzip.open(filename)
	D=cPickle.load(file)
	file.close()
	return D
    def getAMI(self, window=200):
	if self._setAMI == False:
	    self.ami = getAMI(self.seq, window)
	    self._setAMI = True
    def getDBP1(self):
	if self._setDBP1 == False:
	    # THIS IS A CHEAP FIX DON'T RELY ON IT
	    refD, refV = self._dbpload("/home/keith/Documents/filesForProgramming/Libraries/kmbio/kmbio/mast1ref.gzip")
	    dbp, wc = build_Fragment_Dictionary_LZ78_av(self.seq, refD, len(refV))
	    self.dbp1 = dbp
	    self._setDBP1 = True
    def getDBP2(self):
	if self._setDBP2 == False:
	    # THIS IS A CHEAP FIX DON'T RELY ON IT
	    refD, refV = self._dbpload("/home/keith/Documents/filesForProgramming/Libraries/kmbio/kmbio/mast2ref.gzip")
	    dbp, wc = build_Fragment_Dictionary_LZ78_av(self.seq, refD, len(refV))
	    self.dbp2 = dbp
	    self._setDBP2 = True
    def getDBP(self, refD, refV):
	dbp, wc = build_Fragment_Dictionary_LZ78_av(self.seq, refD, len(refV))
	self.dbp = dbp
    def getDict(self):
	if self._setDict == False:
	    self.dict = getDict(self.seq)
	    self._setDict = True
    def getKmer(self, k=0, kobject=None):
	if self._setKmer == False and k==0:
	    k = 7
	if (kobject == None) and (k!=0):
	    kobject=dnakmer(k)
	# Because I have a default kmer, I need to figure out when
	# No k is provided and I should use the default,
	# or when I'm just checking ot see if a k was usedlaksjdflkajsdf
	if (self._setKmer == False) or ((self._kVal != k) and (k != 0)):
	    self.kmer = kobject.profile(self.seq)
	    self._setKmer = True
	    self._kVal = k
	    
    def profileIt(self):
	self.getDBP1()
	self.getDBP2()
	self.getAMI()
	self.getDict()
	self.getKmer()

class setSequences(object):
    def __init__(self, seqIDs, Seqs):
	self.seqs = []
	self.seqIDs = seqIDs
	self._setDM = False
	for i in range(len(seqIDs)):
	    seqIDs[i] = re.sub(r",().:;'\t", "",seqIDs[i])
	    self.seqs.append(sequence(seqIDs[i], Seqs[i]))
    def getAMI(self, window=200):
	for i in range(len( self.seqs)):
	    self.seqs[i].getAMI(window)
    def getDBP1(self):
	for i in range(len( self.seqs)):
	    self.seqs[i].getDBP1()
    def getDBP2(self):
	for i in range(len( self.seqs)):
	    self.seqs[i].getDBP2()
    def getDBP(self, refD, refV):
	for i in range(len( self.seqs)):
	    self.seqs[i].getDBP(refD, refV)
    def getKmer(self, k=0):
	kobject = dnakmer(k)
	for i in range(len(self.seqs)):
	    self.seqs[i].getKmer(k, kobject)
    def getDict(self):
	for i in range(len( self.seqs)):
	    self.seqs[i].getDict()
    def getProfiles(self):
	if len(self.seqs) > 20:
	    print "This will take some time.."
	for i in range(len(self.seqs)):
	    print i
	    self.seqs[i].profileIt()
    def distMat(self, profile, distType):
	dm = [[0.0 for y in range(len(self.seqs))] for x in range(len(self.seqs))]
	labels = []
	if profile.upper() == "DBP1":
	    if distType.upper() == "C":
		for i in range(len( self.seqs)):
		    self.seqs[i].getDBP1()
		    labels.append(self.seqs[i].header)
		    for j in range(i):
			if j != i:
			    dm[i][j] = correlation(self.seqs[i].dbp1, self.seqs[j].dbp1)
			    dm[j][i] = dm[i][j]
	    elif distType.upper() == "E":
		for i in range(len( self.seqs)):
		    self.seqs[i].getDBP1()
		    labels.append(self.seqs[i].header)
		    for j in range(i):
			if j != i:
			    dm[i][j] = euclidean(self.seqs[i].dbp1, self.seqs[j].dbp1)
			    dm[j][i] = dm[i][j]
	elif profile.upper() == "DBP2":
	    if distType.upper() == "C":
		for i in range(len( self.seqs)):
		    self.seqs[i].getDBP2()
		    labels.append(self.seqs[i].header)
		    for j in range(i):
			if j != i:
			    dm[i][j] = correlation(self.seqs[i].dbp2, self.seqs[j].dbp2)
			    dm[j][i] = dm[i][j]
	    elif distType.upper() == "E":
		for i in range(len( self.seqs)):
		    self.seqs[i].getDBP2()
		    labels.append(self.seqs[i].header)
		    for j in range(i):
			if j != i:
			    dm[i][j] = euclidean(self.seqs[i].dbp2, self.seqs[j].dbp2)
			    dm[j][i] = dm[i][j]
	elif profile.upper() == "DBP":
	    if distType.upper() == "C":
		for i in range(len( self.seqs)):
		    labels.append(self.seqs[i].header)
		    for j in range(i):
			if j != i:
			    dm[i][j] = correlation(self.seqs[i].dbp, self.seqs[j].dbp)
			    dm[j][i] = dm[i][j]
	    elif distType.upper() == "E":
		for i in range(len( self.seqs)):
		    labels.append(self.seqs[i].header)
		    for j in range(i):
			if j != i:
			    dm[i][j] = euclidean(self.seqs[i].dbp, self.seqs[j].dbp)
			    dm[j][i] = dm[i][j]
	elif profile.upper() == "AMI":
	    if distType.upper() == "C":
		for i in range(len( self.seqs)):
		    self.seqs[i].getAMI()
		    labels.append(self.seqs[i].header)
		    for j in range(i):
			if j != i:
			    dm[i][j] = correlation(self.seqs[i].ami, self.seqs[j].ami)
			    dm[j][i] = dm[i][j]
	    elif distType.upper() == "E":
		for i in range(len( self.seqs)):
		    self.seqs[i].getAMI()
		    labels.append(self.seqs[i].header)
		    for j in range(i):
			if j != i:
			    dm[i][j] = euclidean(self.seqs[i].ami, self.seqs[j].ami)
			    dm[j][i] = dm[i][j]
	elif profile.upper() == "KMER":
	    if distType.upper() == "C":
		for i in range(len( self.seqs)):
		    self.seqs[i].getKmer()
		    labels.append(self.seqs[i].header)
		    for j in range(i):
			if j != i:
			    dm[i][j] = correlation(self.seqs[i].kmer, self.seqs[j].kmer)
			    dm[j][i] = dm[i][j]
	    elif distType.upper() == "E":
		for i in range(len( self.seqs)):
		    self.seqs[i].getKmer()
		    labels.append(self.seqs[i].header)
		    for j in range(i):
			if j != i:
			    dm[i][j] = euclidean(self.seqs[i].kmer, self.seqs[j].kmer)
			    dm[j][i] = dm[i][j]
	elif profile.upper() == "DICT":
	    if distType.upper() == "C":
		for i in range(len( self.seqs)):
		    self.seqs[i].getDict()
		    labels.append(self.seqs[i].header)
		    for j in range(i):
			if j != i:
			    dm[i][j] = basicGramD(self.seqs[i].dict, self.seqs[j].dict)
			    dm[j][i] = dm[i][j]
	    elif distType.upper() == "E":
		for i in range(len( self.seqs)):
		    self.seqs[i].getDict()
		    labels.append(self.seqs[i].header)
		    for j in range(i):
			if j != i:
			    dm[i][j] = basicGramD(self.seqs[i].dict, self.seqs[j].dict)
			    dm[j][i] = dm[i][j]
	else:
	    raise SyntaxError("Bad input to distMat")
	    exit()
	print "DONE"
	self.dm = dm
	self._setDM = True
	return #dm, labels


class dnakmer(object):
    def __init__(self, size):
	self.vect = [0.0]*(4**size)
	self.k = size
	self.map = self._mapping(size)
    def _mapping(self, k):
	mapStr = 'A'*k
	index = 0
	dnaSet = ['A', 'C', 'G', 'T']
	controlFlow = [0]*k
	strList = ['A']*k
	Dmap = {}
	for i in range(4**k):
	    Dmap[''.join(strList)] = i
	    controlFlow[-1] += 1
	    for j, v in reversed(list(enumerate(controlFlow))):
		if controlFlow[j] == 4:
		    controlFlow[j] =0
		    controlFlow[j-1] += 1
		strList[j] = dnaSet[controlFlow[j]]
	return Dmap
    def profile(self, seq):
	seq = seq.upper()
	vect = [0.0]*4**self.k
	for i in range(len(seq)-self.k+1):
	    if seq[i:i+self.k] in self.map:
		vect[self.map[seq[i:i+self.k]]] += 1.
	vect = numpy.array(vect)/float(len(seq)-self.k+1)
	return vect
	

def getTreeUPGMASet(sampleSet):
    def formatDM(dmraw, labels):
	matrix = []
	for i in range(len(dmraw)):
	    subM = []
	    for j in range(0, i+1):
		subM.append(dmraw[i][j])
	    matrix.append(subM)
	#print matrix
	return _DistanceMatrix(labels, matrix)
    # Take in DM, and Labels
    # format proper for biopy
    dmraw = sampleSet.dm
    labels = sampleSet.seqIDs
    # Strip commas 
    b=",():;'\t"
    for i in range(len(labels)):
	labels[i] = ''.join(c for c in labels[i] if c not in b)
    dm = formatDM(dmraw, labels)
    # Run biopy
    # take tree file out
    constructor = DistanceTreeConstructor()
    treeUPGMA = constructor.upgma(dm)
    Phylo.write(treeUPGMA, 'temp.nwk', 'newick')
    #print treeNJ
    # pass to tree generator: ete3
    tree = Tree('temp.nwk', format=1)
    tree.show()
    return

def getTreeNJSet(sampleSet):
    def formatDM(dmraw, labels):
	matrix = []
	for i in range(len(dmraw)):
	    subM = []
	    for j in range(0, i+1):
		subM.append(dmraw[i][j])
	    matrix.append(subM)
	#print matrix
	return _DistanceMatrix(labels, matrix)
    # Take in DM, and Labels
    # format proper for biopy
    dmraw = sampleSet.dm
    labels = sampleSet.seqIDs
    # Strip commas 
    b=",().:;'\t"
    for i in range(len(labels)):
	labels[i] = ''.join(c for c in labels[i] if c not in b)
    dm = formatDM(dmraw, labels)
    # Run biopy
    # take tree file out
    constructor = DistanceTreeConstructor()
    treeNJ = constructor.nj(dm)
    Phylo.write(treeNJ, 'temp.nwk', 'newick')
    #print treeNJ
    # pass to tree generator: ete3
    tree = Tree('temp.nwk', format=1)
    tree.show()
    return

def nSpectSet(sampleSet):
    dm = sampleSet.dm
    labels = sampleSet.seqIDs
    count = len(labels)
    print 1
    dmf = open('tempDM.dm', 'w')
    lbf = open('tempLabel.txt', 'w')
    dspf = open('tempDisplay.txt', 'w')
    dmf.write(str(count)+"\n")
    for i in range(count):
	lbf.write(str(labels[i])+"\n")
	dspf.write("1 0 4\n")
	dmf.write("seq" + str(i) + "\t")
	for j in range(len(dm[i])):
	    #print dmf[i][j]
	    dmf.write(str.format("{0:.10f}",dm[i][j]))
	    if (j == count-1):
		dmf.write('\n')
	    else:
		dmf.write('\t')
    dmf.close()
    lbf.close()
    dspf.close()
    print 2
    '''
    http://stackoverflow.com/questions/5137726/creating-permanent-executable-aliases
    Ok basically i'm going to make nSpect a sys wide program
    ...Whoop...
    Nevermind..
    This then!
    http://stackoverflow.com/questions/431684/how-do-i-cd-in-python
    One day...
    http://stackoverflow.com/questions/12590262/how-can-i-get-python-subprocess-stdout-and-stderr-to-show-up-in-visual-studios
    p=subprocess.Popen("./nSpect -i "+oldPath+"/tempDM.dm -l "+oldPath+"/tempLabel.txt -d "+oldPath+"/tempDisplay.txt", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout, stderr) = 	p.communicate()
    print stdout
    print stderr
    '''
    cdStr = ""
    oldPath = os.path.dirname(os.path.abspath('tempDM.dm')) # reset tempDM.dm to a varible
    # This is so ugly and hackish. But hey, it works. :D
    os.chdir("/home/keith/Documents/filesForProgramming/Libraries/kmbio/kmbio/")
    #os.system("./nSpect -i "+oldPath+"/tempDM.dm -l "+oldPath+"/tempLabel.txt -d "+oldPath+"/tempDisplay.txt")
    p = Popen("./nSpect -i "+oldPath+"/tempDM.dm -l "+oldPath+"/tempLabel.txt -d "+oldPath+"/tempDisplay.txt", stdout = PIPE, stderr = STDOUT, shell = True)
    for line in p.stdout:
	print line.replace('\n','')
    while True:
	try:
	    print p.stdout.next().replace('\n', '')
	except StopIteration:
	    break
    os.chdir(oldPath)


def distMat(profiles, distType="C"):
    dm = [[0.0 for y in range(len(profiles))] for x in range(len(profiles))]
    if distType.upper() == "C":
	for i in range(len(profiles)):
	    for j in range(i):
		if j != i:
		    dm[i][j] = correlation(profiles[i].T, profiles[j].T)
		    dm[j][i] = dm[i][j]
    elif distType.upper() == "E":
	for i in range(len(profiles)):
	    for j in range(i):
		if j != i:
		    dm[i][j] = euclidean(profiles[i], profiles[j])
		    dm[j][i] = dm[i][j]
    #print dm
    return dm

def nSpect(dm, labels=[], display=[], labelFileName="", displayFileName=""):
    labelVals = False
    displayVals = False
    lfn = False
    dfn = False
    dmFilePath = False
    if not labels==[]:
	labelVals = True
    if not display==[]:
	displayVals = True
    if not labelFileName =="":
	lfn = True
    if not displayFileName =="":
	dfn = True
    if type(dm)==str:
	dmFilePath = True
    
    #print dm[0][2]
    if dmFilePath == False:
	# Assume you're handed a matrix if no fp is given (ie type == str)
	count = len(dm)
	dmf = open('nSpectDM.dm', 'w')
	dmf.write(str(count)+"\n")
	for i in range(count):
	    dmf.write("seq" + str(i) + "\t")
	    for j in range(len(dm[i])):
		#print dmf[i][j]
		dmf.write(str.format("{0:.10f}",dm[i][j]))
		if (j == count-1):
		    dmf.write('\n')
		else:
		    dmf.write('\t')
	dmf.close()

    if labelVals == True:
	lbf = open('nSpectLabel.txt', 'w')
	for i in range(len(labels)):
	    lbf.write(str(labels[i])+"\n")
	lbf.close()


    if displayVals == True:
	dspf = open('nSpectDisplay.txt', 'w')
	for i in range(len(display)):
	    dspf.write(str(display[i][0])+" " +str(display[i][1])+" " +str(display[i][2])+ "\n" )
	dspf.close()




    '''
    http://stackoverflow.com/questions/5137726/creating-permanent-executable-aliases
    Ok basically i'm going to make nSpect a sys wide program
    ...Whoop...
    Nevermind..
    This then!
    http://stackoverflow.com/questions/431684/how-do-i-cd-in-python
    One day...
    http://stackoverflow.com/questions/12590262/how-can-i-get-python-subprocess-stdout-and-stderr-to-show-up-in-visual-studios
    p=subprocess.Popen("./nSpect -i "+oldPath+"/tempDM.dm -l "+oldPath+"/tempLabel.txt -d "+oldPath+"/tempDisplay.txt", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout, stderr) = 	p.communicate()
    print stdout
    print stderr
    '''
    cmdStr = ""
    oldPath = os.getcwd()# os.path.dirname(os.path.abspath(os.getcwd())) # reset tempDM.dm to a varible
    
    if dmFilePath:
	dmFileNamePath = os.path.dirname(os.path.abspath(dm))
	cmdStr = cmdStr + "-i " + dmFileNamePath+ "/" +dm + " "
    else:
	cmdStr = cmdStr + "-i " + oldPath + "/nSpectDM.dm "
    if labelVals:
	cmdStr = cmdStr + "-l " + oldPath + "/nSpectLabel.txt "
    elif lfn:
	labelFileNamePath = os.path.dirname(os.path.abspath(labelFileName))
	cmdStr = cmdStr + "-l " + labelFileNamePath + "/" + labelFileName+ " "
    if displayVals:
	cmdStr = cmdStr + "-d " + oldPath + "/nSpectDisplay.txt "
    elif lfn:
	displayFileNamePath = os.path.dirname(os.path.abspath(displayFileName))
	cmdStr = cmdStr + "-d " + displayFileNamePath + "/" + displayFileName + " "
	


    # This is so ugly and hackish. But hey, it works. :D
    os.chdir("/home/keith/Documents/filesForProgramming/Libraries/kmbio/kmbio/")
    #os.system("./nSpect -i "+oldPath+"/tempDM.dm -l "+oldPath+"/tempLabel.txt -d "+oldPath+"/tempDisplay.txt")
    #print "./nSpect " + cmdStr
    p = Popen("./nSpect " + cmdStr, stdout = PIPE, stderr = STDOUT, shell = True)
    for line in p.stdout:
	print line.replace('\n','')
    while True:
	try:
	    print p.stdout.next().replace('\n', '')
	except StopIteration:
	    break
    os.chdir(oldPath)

def upgmaTree(dmraw, labels):
    def formatDM(dmraw, labels):
	matrix = []
	for i in range(len(dmraw)):
	    subM = []
	    for j in range(0, i+1):
		subM.append(dmraw[i][j])
	    matrix.append(subM)
	#print matrix
	return _DistanceMatrix(labels, matrix)
    # Take in DM, and Labels
    # format proper for biopy
    # Strip commas 
    b=",().:;'\t"
    for i in range(len(labels)):
	labels[i] = ''.join(c for c in labels[i] if c not in b)
    dm = formatDM(dmraw, labels)
    # Run biopy
    # take tree file out
    constructor = DistanceTreeConstructor()
    treeUPGMA = constructor.upgma(dm)
    Phylo.write(treeUPGMA, 'temp.nwk', 'newick')
    #print treeNJ
    # pass to tree generator: ete3
    tree = Tree('temp.nwk', format=1)
    tree.show()
    return

def njTree(dmraw, labels):
    def formatDM(dmraw, labels):
	matrix = []
	for i in range(len(dmraw)):
	    subM = []
	    for j in range(0, i+1):
		subM.append(dmraw[i][j])
	    matrix.append(subM)
	#print matrix
	return _DistanceMatrix(labels, matrix)
    # Take in DM, and Labels
    # format proper for biopy
    # Strip commas 
    b=",().:;'\t"
    for i in range(len(labels)):
	labels[i] = ''.join(c for c in labels[i] if c not in b)
    dm = formatDM(dmraw, labels)
    # Run biopy
    # take tree file out
    constructor = DistanceTreeConstructor()
    treeNJ = constructor.nj(dm)
    Phylo.write(treeNJ, 'temp.nwk', 'newick')
    #print treeNJ
    # pass to tree generator: ete3
    tree = Tree('temp.nwk', format=1)
    tree.show()
    return
def labelLoad(filename):
    fl = open(filename, 'r')
    labels = []
    for line in fl:
	if len(line) > 0:
	    labels.append(line.strip())
    return labels

def loadDM_nSpect(filename):
    # This file assumes nSpect Valid DM
    fl = open(filename, 'r')
    dm = []
    for line in fl:
	line = line.split('\t')
	if len(line) > 1:
	    temp = line[1:]
	    for i in range(len(temp)):
		temp[i] = float(temp[i])
	    dm.append(temp)
    return numpy.array(dm)
def loadDM_raw(filename):
    # This file assumes nSpect Valid DM
    fl = open(filename, 'r')
    dm = []
    for line in fl:
	line = line.strip().split('\t')
	if len(line) > 1:
	    temp = line
	    for i in range(len(temp)):
		temp[i] = float(temp[i])
	    dm.append(temp)
    return numpy.array(dm)

def correlation(a, b):
    return 1 - math.fabs(scipy.stats.pearsonr(a, b)[0])
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
    if denom == 0:
	denom = .000000001
    distance = float(dotProduct)/denom
    correlationDistance = 1 - distance 
    return correlationDistance
    '''

def euclidean(a, b) :
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
    '''
    if ( euclideanDistance != 0) :
        euclideanDistance = 1 / euclideanDistance
    else :
        euclideanDistance = 1
    '''
    return euclideanDistance
def basicGramD(D1, D2):
    union = 0
    intersect = 0
    for key in D1:
	if key in D2:
	    intersect += 1
	    #D2.pop(key, None) # Adds = opperations it removes
	union += 1
    for key in D2:
	if key not in D1:
	    union += 1
    
    '''
    Time tests:
    set opperations:	delt t =05:20
	Getting Distance  2016-03-22 14:19:10
	Done  2016-03-22 14:24:30
    for loop w/o pop:	delt t =02:26
	Getting Distance  2016-03-22 14:26:16
	Done  2016-03-22 14:28:42
    for loop w/ pop:	delt t =02:25
	Getting Distance  2016-03-22 14:31:15
	Done  2016-03-22 14:33:40

    '''


    #num = len(set(D1.keys())&set(D2.keys()))
    #denom = len(set(D1.keys())|set(D2.keys()))
    return 1. - intersect/float(union)

def dnaDataConvert(name):
#Converts string of letters (ATCG) into an
#array of numbers and returns it as "dataConverted"
        # Base A becomes 0
        # Base T becomes 1
        # Base G becomes 2
        # Base C becomes 3
    name = name.replace('A', '0')
    name = name.replace('a', '0')
    name = name.replace('T', '1')
    name = name.replace('t', '1')
    name = name.replace('G', '2')
    name = name.replace('g', '2')
    name = name.replace('C', '3')
    name = name.replace('c', '3')
    name = name.replace('N', '4')
    name = name.replace('R', '4')
    name = name.replace('K', '4')    
    name = name.replace('Y', '4')
    name = name.replace('W', '4')
    name = name.replace('M', '4')
    name = name.replace('S', '4')    
    name = name.replace('B', '4')
    name = name.replace('V', '4')
    name = name.replace('D', '4')
    name = name.replace('E', '4')
    name = name.replace('F', '4')
    name = name.replace('H', '4')
    name = name.replace('I', '4')    
    name = name.replace('J', '4')
    name = name.replace('L', '4')
    name = name.replace('O', '4')
    name = name.replace('P', '4')  
    name = name.replace('Q', '4')
    name = name.replace('U', '4')
    name = name.replace('X', '4')
    name = name.replace('Z', '4')
    name = name.replace('b', '4')
    name = name.replace('d', '4')
    name = name.replace('e', '4')
    name = name.replace('f', '4')
    name = name.replace('h', '4')
    name = name.replace('i', '4')
    name = name.replace('j', '4')
    name = name.replace('k', '4')
    name = name.replace('l', '4')
    name = name.replace('m', '4')
    name = name.replace('n', '4')
    name = name.replace('o', '4')
    name = name.replace('p', '4')
    name = name.replace('q', '4')
    name = name.replace('r', '4')
    name = name.replace('s', '4')
    name = name.replace('u', '4')
    name = name.replace('v', '4')
    name = name.replace('w', '4')
    name = name.replace('x', '4')
    name = name.replace('y', '4')
    name = name.replace('z', '4')
    baseList = list(name)
    baseList = [ int(a) for a in baseList ] 
    length = len(baseList)        
    dataConverted = baseList
    return dataConverted
def amiMatrixProbability(dataConverted) :
    fullBaseList = dataConverted
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
        if ((length != 0) and (a != 4)):
            a = float(a) / length
        independentProb[0][i] = a
        a = 0
        i += 1
    return independentProb
def amiJointProbability(dataConverted, window) :
    fullBaseList = dataConverted
    length = len(fullBaseList)
    independentProb = amiMatrixProbability(dataConverted)
    #window
    k = 0
    amiProfileList = []
    #window = 64

    for k in range(1, int(window) + 2) :
        i = 0
        j = 0
        jointProb =  [ [ 0. for a in range(5) ] for b in range(5) ]
        
        for window in range(length - k) :
            i = fullBaseList[window]
            j = fullBaseList[window + k]
            jointProb[i][j] += 1

        i = 0
        j = 0
        for i in range(4) :
            for j in range(4) :
                if ((length-k) != 0) :
                    jointProb[i][j] = (jointProb[i][j]) / (length - k) 

        #Log Math time
        iOfK = 0
        logNum = 0
        for i in range(4) :
            for j in range(4) :
                if ((independentProb[0][i] != 0)and(independentProb[0][j] !=0)) : 
                    logNum = (float(jointProb[i][j])/ (float(independentProb[0][i]) * independentProb[0][j]) )                
                if (logNum != 0) :
                    iOfK += (float(jointProb[i][j]) * math.log(logNum))
                
        amiProfileList.append(iOfK)
    #*****#End for loop
    return amiProfileList

def getAMI(seq, window) :   
    dnaConvertedData = dnaDataConvert(seq)
    amiProfile = amiJointProbability(dnaConvertedData, window)
    return amiProfile

def getDict(s): # LZ78 algorithm looking for words at each position
    s=s.upper()     # make the sequence uppercase
    D={'A':0,'C':0,'G':0,'T':0}     # initial dictionary

    for i in xrange(len(s)):        # do for the whole sequence
        word_len=1                  # start with 1 letter word
        while s[i:i+word_len] in D: #if we have this word in the dictionary
	    #if s[i+word_len-1:i+word_len] == "N":
		#break
            D[s[i:i+word_len]] +=1  # counting how many times we look for this word in the dictionary
            word_len +=1            # we had the word in the dictionary, now look for a longer one with an extra letter at the end
            if (i+word_len)>=len(s):    #if we reached the end stop the process
                break
        if (i+word_len)<=len(s):  # if we are at the end of the sequence, don't to anything
	    if s[i+word_len-1:i+word_len] !="N": # Make sure unknown characters aren't factored in
                D[s[i:i+word_len]] = 1   # this is a new word, add it to the dictionarym and say it we encountered it once
    
    return D    # report the dictionary to the program
# *******************

def ref_Gram_Vect(s): # Modified LZ78 algorithm looking for words at each position
    s=s.upper()     # make the sequence uppercase
    D={'A':0,'C':1,'G':2,'T':3}     # initial dictionary, value point to vector location
    vector = [0, 0, 0, 0, 0]
    wordLoc = 0
    wordCount = 5
    refD = {0:'A', 1:'C', 2:'G', 3:'T'} 
    

    errorCount = 0
    orderedWords = ['A','C','G','T']
    stepCount = 1
    i = 0
    while i <= len(s):
    #for i in xrange(len(s)):        # do for the whole sequence
        word_len=1                  # start with 1 letter word
        while s[i:i+word_len] in D: #if we have this word in the dictionary            
            wordLoc = D[s[i:i+word_len]]  # reassigning 
	    #print "       " + str(stepCount) + ' & "' + s[i:i+word_len] + '" & ['+str(i)+':'+str(i+word_len-1)+'] & \{'+str(', '.join(orderedWords))+'\} & Extend Window \\\\ \n       \hline}'						
	    vector[wordLoc] += 1  #incrementing vectors value
            word_len +=1            # we had the word in the dictionary, now look for a longer one with an extra letter at the end
	    stepCount += 1
            if (i+word_len)>=len(s):    #if we reached the end stop the process
                break

        if (i+word_len)<=len(s):  # if we are at the end of the sequence, don't to anything
	    #print "       " + str(stepCount) + ' & "' + s[i:i+word_len] + '" & ['+str(i)+':'+str(i+word_len-1)+'] & \{'+(', '.join(orderedWords))+'\} & Add Term to Dictionary \\\\ \n       \hline}'
	    stepCount += 1
	    if s[i+word_len-1:i+word_len] !="N": # Make sure unknown characters aren't factored in
                D[s[i:i+word_len]] = wordCount   # this is a new word, add it to the dictionarym and say it we encountered it once

            #D[s[i:i+word_len]] = wordCount   # this is a new word, add it to the dictionary and say it we encountered it once
	    #refD[wordCount] = s[i:i+word_len]
	    wordCount += 1
	    vector.append(1)
	    orderedWords.append(s[i:i+word_len])
	    #print "       " + str(stepCount) + ' & "' + s[i:i+word_len] + '" & ['+str(i)+':'+str(i+word_len-1)+'] & \{'+(', '.join(orderedWords))+'\} & Reset Term \\\\ \n       \hline}'						
	    stepCount += 1
	    i = i+word_len-2
	    if ( (len(D.keys()) != len(vector)-1) and (errorCount < 2) ):
		print '***********'
		print "ERROR OCCURED"
		print wordCount
		print len(D.keys()), len(vector)
		print D[s[i:i+word_len]]
		print s[i:i+word_len]
		print len(s)
		print i+word_len
		print '***********'
		errorCount += 1
	i += 1
        #if (i+word_len)==len(s): 
	    
    #print D
    return D, vector, refD    # report the dictionary to the program

def gram_Vect(s, D, vecLen): # Modified LZ78 algorithm looking for words at each position
    s=s.upper()     # make the sequence uppercase
    #D={'A':0,'C':1,'G':2,'T':3, 'N':4}     # initial dictionary, value point to vector location
    vector = [0] * vecLen
    wordLoc = 0
    wordCount = 4
    #print D
    #print type(s[0]), D['A']
    i = 0
    while i <= len(s):        # do for the whole sequence
        word_len=1                  # start with 1 letter word
        while s[i:i+word_len] in D: #if we have this word in the dictionary            
            wordLoc = D[s[i:i+word_len]]  # reassigning 
	    vector[wordLoc] += 1  #incrementing vectors value
            word_len +=1            # we had the word in the dictionary, now look for a longer one with an extra letter at the end
	    wordCount += 1
	    #print "YOUOU"
            if (i+word_len)>=len(s):    #if we reached the end stop the process
                break
	if (i+word_len)>=len(s):    #if we reached the end stop the process
            break
	#print "YOUOU"
	i = i+word_len
    for i in range(len(vector)):
	vector[i] = float(vector[i]/float(len(s)) )
    #print vector
    return vector, wordCount    # report the fragVector to the program

# *******************
def build_Master_Dictionary_LZ78_av(s): # Modified LZ78 algorithm looking for words at each position
    s=s.upper()     # make the sequence uppercase
    D={'A':0,'C':1,'G':2,'T':3}     # initial dictionary, value point to vector location
    vector = [0, 0, 0, 0, 0]
    wordLoc = 0
    wordCount = 5
    refD = {0:'A', 1:'C', 2:'G', 3:'T', 4:'N'} 

    errorCount = 0
    for i in xrange(len(s)):        # do for the whole sequence
        word_len=1                  # start with 1 letter word
        while s[i:i+word_len] in D: #if we have this word in the dictionary            
            wordLoc = D[s[i:i+word_len]]  # reassigning 
	    vector[wordLoc] += 1  #incrementing vectors value
            word_len +=1            # we had the word in the dictionary, now look for a longer one with an extra letter at the end
            if (i+word_len)>=len(s):    #if we reached the end stop the process
                break
        if (i+word_len)<len(s):  # if we are at the end of the sequence, don't to anything
            D[s[i:i+word_len]] = wordCount   # this is a new word, add it to the dictionary and say it we encountered it once
	    refD[wordCount] = s[i:i+word_len]
	    wordCount += 1
	    vector.append(1)
	    if ( (len(D.keys()) != len(vector)) and (errorCount < 2) ):
		print '***********'
		print "ERROR OCCURED"
		print wordCount
		print D[s[i:i+word_len]]
		print s[i:i+word_len]
		print len(s)
		print i+word_len
		print '***********'
		errorCount += 1
    return D, vector, refD    # report the dictionary to the program

def build_Fragment_Dictionary_LZ78_av(s, D, vecLen): # Modified LZ78 algorithm looking for words at each position
    s=s.upper()     # make the sequence uppercase
    #D={'A':0,'C':1,'G':2,'T':3, 'N':4}     # initial dictionary, value point to vector location
    vector = [0] * vecLen
    wordLoc = 0
    wordCount = 4

    for i in xrange(len(s)):        # do for the whole sequence
        word_len=1                  # start with 1 letter word
        while (s[i:i+word_len] in D) and s[i+word_len-1] != "N": #if we have this word in the dictionary            
            wordLoc = D[s[i:i+word_len]]  # reassigning 
	    vector[wordLoc] += 1  #incrementing vectors value
            word_len +=1            # we had the word in the dictionary, now look for a longer one with an extra letter at the end
	    wordCount += 1
            if (i+word_len)==len(s)+1:    #if we reached the end stop the process
                break
	if (i+word_len)==len(s):    #if we reached the end stop the process
            break
    for i in range(len(vector)):
	vector[i] = float(vector[i]/float(len(s)) )
    vector[4] = 0.
    
    return vector, wordCount    # report the fragVector to the program


def fna_read(name):
    fid = open(name)
    header = fid.readline()
    dna = fid.read().replace("\n","")
    dna = dna.replace("\r","")
    fid.close()
    return header, dna


def multifna_read(name):
    fragList = []
    headerList = []
    fid = open(name)
    a = fid.readline()
    headerList.append(a.strip())
    s = ''
    for line in fid:
        if line[0]== '>':
	    headerList.append(line.strip())
            fragList.append(s.upper())
            s=''
        else:
	    temp = line.strip()
	    temp = temp.upper()
	    for i in range(len(temp)):
		if temp[i] not in ["A", "T", "G", "C", "N"]:
		    temp = temp.split(temp[i])
		    temp = "N".join(temp)
            s = s+temp
    fid.close()
    fragList.append(s)
    #print len(fragList)
    return headerList, fragList
def read_Fasta(name) :
    #Function reads data from file passed in (name)
        #returns an array with the dna as Strings,
        #and an array with the names of each section
    seqS =[]
    seqName = []
    #Ex: BioPython SeqIO
    handle = open(name, "rU")
    for record in SeqIO.parse(handle, "fasta") :
        seqS.append(str(record.seq))
        seqName.append(str(record.description))
    handle.close()
    return seqName, seqS

