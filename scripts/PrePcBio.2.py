#!/usr/bin/env python3
import os
import sys
#This is to make cBioPortal.db.1.csv into unique variants file with mutiple tumor types...
#V2:Is to re-align the tumor type more accurataley using a token of group-pmid-tumor_site(universal)... 
try:
	if len(sys.argv) <3:
		print ('USAGE:python cBioPortal.db.1.csv Tumor_Type.*.csv')
		exit(0)
except:
	pass
#WHSC1:ENST00000382891.5:c.3295G>A
fHts = open (sys.argv[2],'r')
lines = fHts.readlines()
types={}
for line in lines:
	tokens  = line.split('\t')
	if tokens[0] =='origin':
		continue
	else:
		parts = tokens[4].split(',')
		for part in parts:
			types[part] = tokens[1]

fH = open (sys.argv[1],'r')
lines = fH.readlines()
D={};l=[]
tt={};pid={};gp={}
header = lines[0];del lines[0]
for line in lines:
	tokens= line.split('\t')
	gene = tokens[0]
	parts = tokens[-5].split(':');l.append(tokens[0]+':'+tokens[-5])
	key = gene +':'+ tokens[-5]#parts[-1]
	tt[key] = []
	del tokens[-1];del tokens[-1];del tokens[-1];tmp = '\t'.join(tokens)
	D[key]= tmp

l = list(set(l))
print ('{}V.{}'.format(len(l),len(D)))
for line in lines:
	tokens= line.split('\t')
	gene = tokens[0]
	parts = tokens[-5].split(':')
	key = gene +':'+ tokens[-5]#parts[-1]
	try:
		tt[key].append(tokens[-3]+':'+tokens[-2]+':'+types[tokens[-1].replace('\n','')])
	except:
		tt[key].append(tokens[-3]+':'+tokens[-2]+':'+tokens[-1].replace('\n',''))
#writing out...
fOut = open ('cBioPortal.db.3.csv','w')
fOut.write('variant\t{}'.format(header))
for var in D:
	tt[var] = list(set(tt[var]))
	fOut.write('{}\t{}\t{}\n'.format(var,D[var],';'.join(tt[var])))	


print ('{}V.{}'.format(len(l),len(D)))






