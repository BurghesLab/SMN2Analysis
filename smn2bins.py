import csv
import subprocess
import re
from collections import defaultdict
import argparse
import sys
import math

parser = argparse.ArgumentParser(description = 'This script will find the copy number of SMN2 over a series of equally sized bins. This script takes a list of sample names as input. There should be one folder for each sample name, with the bam files inside. By default, bam should be named "samplename".star.dedup.bam, though this can be change with the --suffix option.')
requiredNamed = parser.add_argument_group('Required arguments')
requiredNamed.add_argument('-s', '--samples', nargs = '+', type = str, help = 'List of sample names, separated spaces', required = True)
requiredNamed.add_argument('-k', '--known', nargs = '+', type = str, help = 'Samples with known copy number. Should be in the form "samplename,copynumber", eg 112,3', required = True)
parser.add_argument('-x', '--suffix', default = '.star.dedup.bam', type = str, help = 'Bam file suffix. ".star.dedup.bam" by default')
parser.add_argument('-d', '--directory', default = './', type = str, help = 'Directory where all folders containing the bam files are found. Each bam should be in its own folder, which should be named the sample name')
parser.add_argument('-b', '--outBins', default = 'Smn2Bins.txt', help = 'Output file name for binned file. If filename has spaces, please use quotations around it. If no file name is specified, output will be saved to Smn2Bins.txt. If the file exists it will be overwriten!')
parser.add_argument('-c', '--outCopyNum', default = 'Smn2CopyNum.txt', help = 'Output file name for copy number file. If filename has spaces, please use quotations around it. If no file name is specified, output will be saved to Smn2CopyNum.txt. If the file exists it will be overwriten!')

args = parser.parse_args(args = None if len(sys.argv) > 1 else ['--help'])

binFile = args.outBins
copyFile = args.outCopyNum
bamFiles = args.samples
directory = args.directory
suffix = args.suffix
knownInput = args.known

known = []
knownCopies = []
for k in knownInput:
	k = k.split(',')
	known.append(k[0])
	knownCopies.append(int(k[1]))


fasta = 'g1kv37_mask_clone35k.fa'
clones = {'1056O6': 150000, '215P15': 35000, '268A15': 174000, '652K3': 178000}

cftrLen = 150000
smn2Len = 35000
plsLen = 140000


cftrBin = 2499
smn2Bin = 580
plsBin = 2360

def sorted_nicely( l ): 

    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

cftr = []
smn2 =[]
pls3 = []
posDict = defaultdict(lambda: defaultdict(dict))
bedList = []
countDict = {}

for bam in bamFiles:
	bamF = directory + bam + '/' + bam + suffix
	for clone in clones:
		clonLen = clones[clone]
		if clone == '268A15':
			binLen = ((clonLen - 32500) / 19 ) - 1
			i = 32500
		else:
			binLen = ((clonLen - 500) / 19 ) - 1
			i = 500

		j = binLen + i


		while j < (clonLen):
			myString = 'samtools view -c ' + bamF + ' ' + clone + ':' + str(i) + '-' + str(j)
			count = subprocess.check_output([myString], stderr = subprocess.STDOUT, shell = True )
			posDict[bam][clone][str(j)] = count
			bedCmd = clone + '\t' + str(i) + '\t' + str(j) + '\t' + clone +':' + str(j) + '\n'
			bedList.append(bedCmd)
			j += binLen
			i += binLen


adjDict = defaultdict(list)

for j, bam in enumerate(known):
	smn = sorted_nicely(posDict[bam]['215P15'].keys())
	cftr = sorted_nicely(posDict[bam]['1056O6'].keys())
	cftr2 = sorted_nicely(posDict[bam]['652K3'].keys())
	pls = sorted_nicely(posDict[bam]['268A15'])
	for i, coord in enumerate(smn):
		smn2Count = posDict[bam]['215P15'][coord]
		cftrCount = posDict[bam]['1056O6'][cftr[i]]
		cftr2Count = posDict[bam]['652K3'][cftr2[i]]
		plsCount = posDict[bam]['268A15'][pls[i]]
		copyNum = (float(smn2Count)/(float(cftrCount) + float(cftr2Count)))/float(knownCopies[j])
		adjDict[bam].append(copyNum)

binDict = defaultdict(list)
adj = []
for bam in known:
	for i, binn in enumerate(adjDict[bam]):
		binDict[i].append(binn)

for binList in binDict.keys():
	factorAve = 0
	for factor in binDict[binList]:
		factorAve += factor
	factorAve = factorAve / float(len(binDict[binList]))
	adj.append(factorAve)





#
#this is for printing out
with open(binFile, 'wt') as output:
	csvOut = csv.writer(output, delimiter = '\t')
	for bam in bamFiles:
		csvOut.writerow([bam])
		smn = sorted_nicely(posDict[bam]['215P15'].keys())
		cftr = sorted_nicely(posDict[bam]['1056O6'].keys())
		cftr2 = sorted_nicely(posDict[bam]['652K3'].keys())
		pls = sorted_nicely(posDict[bam]['268A15'])
		for i, coord in enumerate(smn):
			smn2Count = posDict[bam]['215P15'][coord]
			cftrCount = posDict[bam]['1056O6'][cftr[i]]
			cftr2Count = posDict[bam]['652K3'][cftr2[i]]
			plsCount = posDict[bam]['268A15'][pls[i]]
			copyNum = (float(smn2Count)/(float(cftrCount) + float(cftr2Count)))/adj[i]
			outLi = [smn2Count.strip(), cftrCount.strip(), cftr2Count.strip(), plsCount.strip(), copyNum]
			csvOut.writerow(outLi)

clones = ['215P15', '1056O6', '652K3', '268A15']

coefficients = {}

maleRatio = []
femaleRatio = []
males = []
#Should give list of known samples on command line, and list of their respective copynum
readsDict = defaultdict(dict)
copyNumDict = defaultdict(list)
sexDict = {}
#Calls samtools to get readcounts of each clone
for bam in bamFiles:
	bamF = directory + bam + '/' + bam + suffix
	for clone in clones:
		myString = 'samtools view -c ' + bamF + ' ' + clone
		count = subprocess.check_output([myString], stderr = subprocess.STDOUT, shell = True )
		readsDict[bam][clone] = count


#Determine the coefficients
#Sex ratio (pls/cftr) usually around 0.28 for males, 0.48 for females
#May vary batch to batch so must determine anew for each run

for bam in readsDict:
	pls = readsDict[bam]['268A15']
	cftr = float(readsDict[bam]['1056O6']) + float(readsDict[bam]['652K3'])
	sexRatio = float(pls) / cftr
	sexDict[bam] = sexRatio

#Find the difference between sexRatio and expected
#The smaller difference determines the sex
#If there is a big deviation it is reported	
	fem = abs(0.59 - sexRatio)
	mal = abs(0.39 - sexRatio)

	if fem > mal:
		male = True
		if mal > 0.1:				
			print 'Warning: '
		else:
			males.append(bam) 
	else:
		male = False
		if fem > 0.1:
			print 'Warning: Error 2'

	if male:
		maleRatio.append(sexRatio)
	else:
		femaleRatio.append(sexRatio)
#Find average then get malefactor, needed when calculating copy number
maleCo = sum(maleRatio)/float(len(maleRatio))
femaleCo = sum(femaleRatio)/float(len(femaleRatio))
maleFactor = maleCo/femaleCo
coefficients['maleFactor'] = maleFactor

#Find coefficients using sample with known copy number

knownCftr = []
knownPls = []
knownAll = []

for i,bam in enumerate(known):
	bamF = directory + bam + '/' + bam + suffix
	
	smn = readsDict[bam]['215P15']
	cftr105 = readsDict[bam]['1056O6']
	cftr65 = readsDict[bam]['652K3']
	pls = readsDict[bam]['268A15']
	cftrCo = (float(smn)/(float(cftr105) + float(cftr65)))/float(knownCopies[i])	

	if bam in males:
		plsCo = (float(smn)/((float(pls)/float(coefficients['maleFactor']))))/float(knownCopies[i])
		allCo = ((float(smn)/((float(pls)/float(coefficients['maleFactor'])) + float(cftr105) + float(cftr65)))/float(knownCopies[i]))
	else:
		plsCo = (float(smn)/float(pls))/float(knownCopies[i])
		allCo = (float(smn)/(float(pls) + float(cftr105) + float(cftr65)))/float(knownCopies[i])

	knownCftr.append(cftrCo)
	knownPls.append(plsCo)
	knownAll.append(allCo)

coefficients['cftr'] = sum(knownCftr)/float(len(knownCftr))
coefficients['pls'] = sum(knownPls)/float(len(knownPls))
coefficients['all'] = sum(knownAll)/float(len(knownAll))

#Find copy number

with open(copyFile, 'wt') as copyoutput:
	csvcopyOut = csv.writer(copyoutput, delimiter = '\t')
	header = ['Sample', '215P15', '1056O6', '652K3', '268A15', 'SexRatio', 'Pls Ratio', 'Cftr Ratio', 'All Ratio', 'SMN Copies']
	csvcopyOut.writerow(header)
	for bam in bamFiles:
		bamF = directory + bam + '/' + bam + suffix

		smn = readsDict[bam]['215P15']
		cftr105 = readsDict[bam]['1056O6']
		cftr65 = readsDict[bam]['652K3']
		pls = readsDict[bam]['268A15']
		sexRatio = sexDict[bam]

		if bam in males:
			plsRatio = (float(smn)/(float(pls)/float(coefficients['maleFactor'])))/float(coefficients['pls'])
			allRatio = ((float(smn)/((float(pls)/float(coefficients['maleFactor'])) + float(cftr105) + float(cftr65)))/float(coefficients['all']))
		else:
			plsRatio = (float(smn)/float(pls))/float(coefficients['pls'])
			allRatio = (float(smn)/(float(pls) + float(cftr105) + float(cftr65)))/float(coefficients['all'])

		cftrRatio = (float(smn)/(float(cftr105) + float(cftr65)))/float(coefficients['cftr'])
		
	
		copyNum = (plsRatio + cftrRatio + allRatio)/float(3)
		outLi = [bam]
		addition = [sexRatio, plsRatio, cftrRatio, allRatio, copyNum]
		decimal = int(str(copyNum % 1)[2])

		if decimal > 6:
			intCopy = int(math.ceil(copyNum))
		if decimal < 7 and decimal > 3:
			print(bam + ' has indeterminate copy number, rounding down')
			intCopy = int(copyNum)
		if decimal < 4:
			intCopy = int(copyNum)

		copyNumDict[intCopy].append(bam)

		for clone in clones:
			outLi.append(readsDict[bam][clone].strip())
			
		outLi = outLi + addition
		csvcopyOut.writerow(outLi)