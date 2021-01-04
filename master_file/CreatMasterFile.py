#!/usr/bin/env python3
import os
import sys
#This is to develop gene census master file based on cancer_gene_census.csv & CCGD_export.db.csv using mapping files TypesMatching.csv & TypesMatching.CCGD.csv...
try:
	if len(sys.argv) < 3:
		print ('USAGE:python cancer_gene_census.csv,CCGD_export.db.csv TypesMatching.csv,TypesMatching.CCGD.csv')
		exit(0)
except:
	pass 
print ('USAGE:python cancer_gene_census.csv,CCGD_export.db.csv TypesMatching.csv,TypesMatching.CCGD.csv')
genesFiles = sys.argv[1].split(',')
typesFiles = sys.argv[2].split(',')
fH1  =open (genesFiles[0],'r')
lines = fH1.readlines()
Genes={}
for line in lines:
	tokens= line.split(',')
	if tokens[0] =='Gene Symbol':
		continue
	else:
		Genes[tokens[0]]={'synonyms':'','type':{'COSMIC':'','CCGD':''},'score':{'all':0, 'aml':0, 'anal':0, 'rectal':0, 'colon':0, 'bladder':0, 'uterine':0, 'bone':0, 'breast':0, 'cns':0, 'esophageal':0, 'gastric':0, 'kaposi':0, 'melanoma':0, 'ovarian':0, 'cervical':0, 'mds':0, 'nhl':0, 'pancreas':0, 'nsclc':0, 'sclc':0, 'hnscc':0}}
		parts = tokens[9].split('+ADs- ')
		part = tokens[4]+'#'+';'.join(parts).replace('+AC0-','-')
		#part = part.replace('+AC0-','-')
		Genes[tokens[0]]['synonyms']=tokens[-1].replace('\n','').replace('+ADs-',';')
		Genes[tokens[0]]['type']['COSMIC']=  part.lower()
print ('NOTICE:Found {} genes in {}...'.format(len(Genes),genesFiles[0]))
fH2 =open (genesFiles[1],'r')
lines = fH2.readlines()
duplicate={'ZMYM2':'ZNF198','SRGAP2':'SRGAP3','TCF4':'TCF7L2','CCL26':'MAX','PATZ1':'ZNF278'}
for line in lines:
	tokens = line.split('\t')
	if tokens[0] =='Gene':
		continue
	else:
		if tokens[0] in Genes:
			tmp = tokens[-2] +'#'+tokens[-1].replace('\n','')			
			Genes[tokens[0]]['type']['CCGD']=tmp
		elif tokens[0] in duplicate:
			tmp = tokens[-2] +'#'+tokens[-1].replace('\n','')
			Genes[duplicate[tokens[0]]]['type']['CCGD']=tmp
		else:
			Genes[tokens[0]]={'synonyms':'','type':{'COSMIC':'','CCGD':[]},'score':{'all':0, 'aml':0, 'anal':0, 'rectal':0, 'colon':0, 'bladder':0, 'uterine':0, 'bone':0, 'breast':0, 'cns':0, 'esophageal':0, 'gastric':0, 'kaposi':0, 'melanoma':0, 'ovarian':0, 'cervical':0, 'mds':0, 'nhl':0, 'pancreas':0, 'nsclc':0, 'sclc':0, 'hnscc':0}}
			tmp = tokens[-2] +'#'+tokens[-1].replace('\n','')			
			Genes[tokens[0]]['type']['CCGD']=tmp
print ('NOTICE:Found {} genes in {} & {}...'.format(len(Genes),genesFiles[0],genesFiles[1]))
print (Genes['PTEN'])
fH3 = open (typesFiles[0],'r')
lines = fH3.readlines()
Ucode=[];COSMIC={}
ignore = ['Not yet','duplicated','NA','???','Sub+AC0-type that not covered']
for line in lines:
	tokens = line.split('\t')
	if tokens[0] =='Census':
		continue
	else:
		Ucode.append(tokens[3])
		if tokens[1] in ignore:
			continue
		else:
			tmp = tokens[0].lower()
			chuncks = tokens[1].split('+ADs-')
			COSMIC[tmp]=chuncks
#			info = Genes[gene]['type']['COSMIC'].split('-')
fH4 = open (typesFiles[1],'r')
lines = fH4.readlines()
CCGD={}
for line in lines:
	tokens= line.split('\t')
	if tokens[0]=='CCGD':
		continue
	elif tokens[0] =='':
		continue
	else:
		if tokens[1] in ignore:
			continue
		else:
			chuncks = tokens[1].split('+ADs-')
			CCGD[tokens[0]]=chuncks
#Dev the score...
for gene in Genes:
	info = Genes[gene]['type']['COSMIC'].split('#')
	if len(info) > 1:
		#if gene == 'MDM4':
		#	print ('^^^^^^^^^^^^^^^^^^^^^^^')
		#	types= info[1].split(';')
		#	print (types)
		#	for typ in types:
		#		if typ in COSMIC:
		#			print (typ)
		#			Genes[gene]['score'][typ]= 1
		#			print (Genes[gene]['score'][typ])
		#	exit(0)
		types= info[1].split(';')#;print (types)
		for typ in types:
			if typ in COSMIC:
				for y in COSMIC[typ]:
					#print ('^^^^^^^^^^');print (Genes[gene]['score'])
					#if gene == 'MDM4':
					#print ('^^^^');print ('{}:{}'.format(y,	typ));print ('^^^^')
					Genes[gene]['score'][y] = 1
					for y in Genes[gene]['score']:
						#print ('^^^^');print ('{}:{}'.format(y,	typ));print ('^^^^')
						if y == typ:
							for u in range (0,len(COSMIC[typ]),++1):
								r=COSMIC[typ][u]#;print (r)
								#print (Genes[gene]['score'][COSMIC[u]])
								Genes[gene]['score'][r] = 1
								#if gene== 'MDM4':
					#			print ('****');print (Genes[gene]['score']);print ('****')
			else:
				continue#print ('WARNING:{}'.format(typ))
			#for t in COSMIC:
			#	for tt in COSMIC[t]:
			#		if typ.lower() ==tt: 
			#			Genes[gene]['score'][tt]= 1
						#print (gene);print(Genes[gene]);exit(0)
	else:#not in COSMIC...
		#alts = Genes[gene]['synonyms'].split(';')
		#print (alts)
		i = Genes[gene]['type']['CCGD'].split('#')
		if len(i) > 1:
			Types = i[1].split(';')
			for Typ in Types:
				if Typ =='':
					continue
				else:
					tparts = Typ.split(':')
					tparts[0] = tparts[0].replace('"','')
					#if gene =='ARNTL':
					#	print ('***********');print(tparts[0])
					#	print(CCGD)
					#	if tparts[0] in CCGD:
					#		print (CCGD[tparts[0]])
					#		for z in CCGD[tparts[0]]:
					#			print (z);Genes[gene]['score'][z]=2
					#			print (Genes[gene]['score'][z])
						#exit(0)
					#for t in CCGD:
						#for ttt in CCGD[t]:
							#print (ttt)
					if tparts[0] in CCGD:
						#print (Typ);print (Types);print(CCGD[t])
						for z in CCGD[tparts[0]]:
							#print (Genes[gene]['score'][z])
							#HERE#print (type(Genes[gene]['score'][z]))
							Genes[gene]['score'][z]=2
						#print (gene);print(Genes[gene]);exit(0)
		else:
			print ('WARNING:Can not find {} information...'.fromat(gene))
			exit(0)
		#print (Genes[gene]);exit(0)
#print('++++++++++')
#print (Genes['PTEN'])
#Writing out...
Ucode=['all','aml','anal','rectal','colon','bladder','uterine','bone','breast','cns','esophageal','gastric','kaposi','melanoma','ovarian','cervical','mds','nhl','pancreas','nsclc','sclc','hnscc']
fOut = open ('CENSUS.MF.db.2.csv','w')
fOut.write('Gene\tSynonmys\tall\taml\tanal\trectal\tcolon\tbladder\tuterine\tbone\tbreast\tcns\tesophageal\tgastric\tkaposi\tmelanoma\tovarian\tcervical\tmds\tnhl\tpancreas\tnsclc\tsclc\thnscc\n')
for gene in Genes:
	fOut.write('{}\t{}\t'.format(gene,Genes[gene]['synonyms']))
	fOut.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(Genes[gene]['score']['all'],Genes[gene]['score']['aml'],Genes[gene]['score']['anal'],Genes[gene]['score']['rectal'],Genes[gene]['score']['colon'],Genes[gene]['score']['bladder'],Genes[gene]['score']['uterine'],Genes[gene]['score']['bone'],Genes[gene]['score']['breast'],Genes[gene]['score']['cns'],Genes[gene]['score']['esophageal'],Genes[gene]['score']['gastric'],Genes[gene]['score']['kaposi'],Genes[gene]['score']['melanoma'],Genes[gene]['score']['ovarian'],Genes[gene]['score']['cervical'],Genes[gene]['score']['mds'],Genes[gene]['score']['nhl'],Genes[gene]['score']['pancreas'],Genes[gene]['score']['nsclc'],Genes[gene]['score']['sclc'],Genes[gene]['score']['hnscc']))
print ('======================')
print (Genes['MDM4'])
exit(0)	
l=[];U=[]
for line in lines:
	tokens = line.split('\t')
	if tokens[0] =='Census':
		continue
	else:
		l.append(tokens[1])
		U.append(tokens[3])
l= list(set(l))
U = list(set(U))
for typ in l:
	tokens = typ.split('+ADs-')
	for token in tokens:
		if token.lower() in U:
			continue
		else:
			print (token)
print (len(l))
#Genes={'COSMIC':'','CCGD':''}

