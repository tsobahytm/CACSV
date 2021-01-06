#!/usr/bin/env python3
import os
import sys
#This is to conert intOgen driver mutations catalog & build db (in CSV)...
#V2:combine mAkeSublists.py & SV-Collection.intOgen.py in a single and more dynamic approach...
#V2.1:to use and associate the universal tumor type to each variant...  

try:
	if len(sys.argv) < 3:
		print ('USAGE:python intOgen_mutation_analysis.csv Tumor_Type.*.csv')
		exit(0)
except:
	pass
	
fHts = open (sys.argv[2],'r')
lines = fHts.readlines()
types={}
for line in lines:
	tokens = line.split('\t')
	if tokens[0] =='origin':
		continue
	else:
		parts = tokens[3].split(',')
		for part in parts:
			types[part]=tokens[1]
fH = open (sys.argv[1],'r')
lines= fH.readlines()
D={}
for line in lines:
	tokens= line.split('\t')
	if tokens[0] == 'sample':
		continue
	else:
		var =tokens[8] +':'+tokens[12]+':'+tokens[13]
		#var =tokens[8] +':'+tokens[13]
	D[var]={'Chr':'','pos':'','ref':'','alt':'','gene':'','region':'','strand':'','transcript':'','gDNA':'','AAmut':'','consequence':'','exon':'','ProteinPosition':'','transcriptAA':'','transcriptExon':'','otherKnownTumor':'','domain':'','generole':'','mutationlocation':'','driverPred':'','TumorType':'','driver':''}
print (len(D))
for line in lines:
	tokens= line.split('\t')
	if tokens[0] =='sample':
		continue
	else:
		var =tokens[8] +':'+tokens[12]+':'+tokens[13]
		#var =tokens[8] +':'+tokens[13]
		if var in D:
			'1	3	4	5	8	10	11	12	14	15	16	17	19	20	21	22	25	28	31	35	36'
			try:
				D[var]['Chr'] += ';'+tokens[1];D[var]['pos'] += ';'+tokens[3];D[var]['ref'] += ';'+tokens[4];D[var]['alt'] += ';'+tokens[5];D[var]['gene'] += ';'+tokens[8];D[var]['region'] += ';'+tokens[10];D[var]['strand'] += ';'+tokens[11];D[var]['transcript'] += ';'+tokens[12];D[var]['gDNA'] += ';'+tokens[14];D[var]['AAmut'] += ';'+tokens[15];D[var]['consequence'] += ';'+tokens[16];D[var]['exon'] += ';'+tokens[17];D[var]['ProteinPosition'] += ';'+tokens[19];D[var]['transcriptAA'] += ';'+tokens[20];D[var]['transcriptExon'] += ';'+tokens[21];D[var]['otherKnownTumor'] += ';'+tokens[22];D[var]['domain'] += ';'+tokens[25];D[var]['generole'] += ';'+tokens[28];D[var]['mutationlocation'] += ';'+tokens[31];D[var]['driverPred'] += ';'+tokens[35];D[var]['TumorType'] += ';'+types[tokens[36]];D[var]['driver'] += ';'+tokens[38].replace('\n','')
			except:
				D[var]['Chr'] += ';'+tokens[1];D[var]['pos'] += ';'+tokens[3];D[var]['ref'] += ';'+tokens[4];D[var]['alt'] += ';'+tokens[5];D[var]['gene'] += ';'+tokens[8];D[var]['region'] += ';'+tokens[10];D[var]['strand'] += ';'+tokens[11];D[var]['transcript'] += ';'+tokens[12];D[var]['gDNA'] += ';'+tokens[14];D[var]['AAmut'] += ';'+tokens[15];D[var]['consequence'] += ';'+tokens[16];D[var]['exon'] += ';'+tokens[17];D[var]['ProteinPosition'] += ';'+tokens[19];D[var]['transcriptAA'] += ';'+tokens[20];D[var]['transcriptExon'] += ';'+tokens[21];D[var]['otherKnownTumor'] += ';'+tokens[22];D[var]['domain'] += ';'+tokens[25];D[var]['generole'] += ';'+tokens[28];D[var]['mutationlocation'] += ';'+tokens[31];D[var]['driverPred'] += ';'+tokens[35];D[var]['TumorType'] += ';'+tokens[36];D[var]['driver'] += ';'+tokens[38].replace('\n','')
			
		else:
			print ('WARNING:Could not find this variant:{}'.format(var))
			exit(0)
Ofile = 'intOgen.db.gene.transcript.cDNA.v2.1.csv'
fOut = open (Ofile,'w')	
fOut.write('gene:cDNA\tChr\tpos\tref\talt\tgene\tregion\tstrand\ttranscript\tgDNA\tAAmut\tconsequence\texon\tProtein Position\ttranscript AA\ttranscript Exon\tother Known Tumor\tdomain\tgene role\tmutation location\tdriver Pred\tTumor Type\tdriver\n')
for var in D:
	fOut.write('{}\t'.format(var))
	for k in D[var]:
		if k == 'driver':
			tokens = D[var][k].split(';')
			tokens = list(set(tokens))
			s = ';'.join(tokens)
			fOut.write('{}\n'.format(s))	
		else:
			tokens = D[var][k].split(';')
			tokens = list(set(tokens))
			s = ';'.join(tokens)
			fOut.write('{}\t'.format(s))
print ('NOTICE:Done writing {} variants into {}...'.format(len(D),Ofile))




