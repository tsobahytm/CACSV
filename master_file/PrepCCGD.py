#!/usr/bin/env python3
import os
import sys
import collections
#This is to prepare CCGD_export.csv to be intergarted into CGC.db.csv...
try:
	if len(sys.argv) <2:
		print ('USAGE:python CCGD_export.csv')
		exit(0)
except:
	pass
def geneKNWL(gene,L):#to get gene-tumor type association that includes number of publications & frequency of predicted effect...
		#testing on example gene...
	#if gene == 'PTEN':
	CGC=[];Pub=[];tType=[]
	for l in L:
		tokens = l.split(',')
		CGC.append(tokens[-8])	
		Pub.append(tokens[-6])
		t= tokens[-2] +':'+ tokens[-6]
		tType.append(t)
	CGC = list(set(CGC))
	if len(CGC) > 1:#to check if a gene has conflict been on COSMIC CGC list... 
		print ('{}---{}'.format(gene,CGC))
	Pub = list(set(Pub))
	counter=collections.Counter(tType)
	fType={}
	for c in counter:
		info = c.split(':')
		fType[info[0]]=[]
	for c in counter:
		info = c.split(':')
		fType[info[0]].append(info[1])
	x=0#to check the quality of frequency calculation...
	for tumor in fType:
		freq = float (len(fType[tumor])) / float (len(Pub)) * 100
		x += freq	
		fType[tumor].append(str(int(freq))+'%')
	if int(x) < 99:
		print ('ERROR:While calculating tumor type frequency per publications...')
		print (x);print (gene);exit(0)

	return CGC[0],Pub,fType
fH = open (sys.argv[1],'r')
lines = fH.readlines()
D={}
for line in lines:
	tokens = line.split(',')
	if tokens[0] == '\n':
		continue
	else:
		if tokens[4] == '"Human Symbol"':
			header = line
		else:
			if tokens[4] =='':#not reported in human yet...
				continue
			else:
				D[tokens[4]]=[]
print (len(D))
for line in lines:
	tokens = line.split(',')
	if tokens[0] == '\n':
		continue
	else:
		if tokens[4] == '"Human Symbol"':
			header = line
		else:
			if tokens[4] =='':#not reported in human yet...
				continue
			else:
				D[tokens[4]].append(line)
Ofile = sys.argv[1].replace('.csv','.db.csv')
fOut = open (Ofile,'w')
fOut.write('Gene\tCGC\tPublications\tCancer type frequency per number of publications\n')
for gene in D:
	CGC, Pub, fType = geneKNWL(gene,D[gene])
	fOut.write('{}\t{}\t{}\t'.format(gene,CGC,len(Pub)))
	for tumor in fType:
		fOut.write('{}:{};'.format(tumor,fType[tumor][-1]))
	fOut.write('\n')

	




	




