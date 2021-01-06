#!/usr/bin/env python
import os
import sys
#This is to get parse allActionableVariants.txt file to converted to db file that will be integrated to our finalDB...
fH = open (sys.argv[1],'r')
lines = fH.readlines()
genes_oncokb={}
for line in lines:
	tokens = line.split('\t')
	if tokens[0] == 'Isoform':
		continue
	else:
		 genes_oncokb[tokens[3]]=[]
for line in lines:
	tokens = line.split('\t')
	if tokens[0] == 'Isoform':
		continue
	else:
		l = [tokens[0],tokens[4],tokens[5],tokens[6],tokens[7],tokens[8]]
		genes_oncokb[tokens[3]].append(l)
print('NOTICE:Found {} genes in {}...'.format(len(genes_oncokb),sys.argv[1]))
fOut = open ('OncoKB.db.1.csv','w')
fOut.write('gene\tENST_ID\talteration\tprotein_change\ttumor_type\tevidence_level\tdrugs\n')
for gene in genes_oncokb:
	#print('{}\n{}'.format(gene,genes_oncokb[gene]))
	for l in genes_oncokb[gene]:
		fOut.write('{}\t{}\n'.format(gene,'\t'.join(l)))
	#exit(0)

#let's investigate tissue types & variants...
#types=[];variants=[]
#for line in lines:
#	tokens = line.split('\t')
#	if tokens[0] == 'Isoform':
#		continue
#	else:
#		variants.append(tokens[4])
#		types.append(tokens[6])
#variants = list(set(variants))
#types = list(set(types))
#l=[]
#for line in lines:
#	tokens = line.split('\t')
	#if 'Exon' in tokens[4]:
	#	l.append(tokens[3])	
	#elif tokens[4] == 'Oncogenic Mutations':
	#	l.append(tokens[3])
#	elif tokens[4] =='Truncating Mutations':
#		l.append(tokens[3])
#print(len(l))
#l = list(set(l));print(len(l))
#print (l)
	
