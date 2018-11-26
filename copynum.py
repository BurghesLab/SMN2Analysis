import csv
import subprocess
import re
from collections import defaultdict
import math

outfile = 'out.txt'
bamFiles = []


clones = ['215P15', '1056O6', '652K3', '268A15']

coefficients = {}
known = []
knownCopies = []
maleRatio = []
femaleRatio = []
males = []

countDict = defaultdict(dict)
copyNumDict = defaultdict(list)
sexDict = {}

for bam in bamFiles:
	bamF = '' + bam + '/' + bam + '' 
	for clone in clones:
		myString = 'samtools view -c ' + bamF + ' ' + clone
		count = subprocess.check_output([myString], stderr = subprocess.STDOUT, shell = True )
		countDict[bam][clone] = count


#Determine the coefficients
#Sex ratio (pls/cftr) usually around 0.28 for males, 0.48 for females
#May vary batch to batch so must determine anew for each run

for bam in countDict:
	pls = countDict[bam]['268A15']
	cftr = float(countDict[bam]['1056O6']) + float(countDict[bam]['652K3'])
	sexRatio = float(pls) / cftr
	print sexRatio
	sexDict[bam] = sexRatio

#Find the difference between sexRatio and expected
#The smaller difference determines the sex
#If there is a big deviation it is reported	
	fem = abs(0.96 - sexRatio)
	mal = abs(0.55 - sexRatio)

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



knownCftr = []
knownPls = []
knownAll = []

for i,bam in enumerate(known):
	bamF = '' + bam + '/' + bam + ''
	
	smn = countDict[bam]['215P15']
	cftr105 = countDict[bam]['1056O6']
	cftr65 = countDict[bam]['652K3']
	pls = countDict[bam]['268A15']
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

with open(outfile, 'wt') as output:
	csvOut = csv.writer(output, delimiter = '\t')
	header = ['Sample', '215P15', '1056O6', '652K3', '268A15', 'SexRatio', 'Pls Ratio', 'Cftr Ratio', 'All Ratio', 'SMN Copies']
	csvOut.writerow(header)
	for bam in bamFiles:
		'' + bam + '/' + bam + '' 

		smn = countDict[bam]['215P15']
		cftr105 = countDict[bam]['1056O6']
		cftr65 = countDict[bam]['652K3']
		pls = countDict[bam]['268A15']
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
			outLi.append(countDict[bam][clone].strip())
			
		outLi = outLi + addition
		csvOut.writerow(outLi)

