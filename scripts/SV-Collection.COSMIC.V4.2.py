#!/usr/bin/env python3
import os
import sys
#This is to create in-house db-COSMIC (in CSV file) ...
#V4:change tumor type to universal per variant as the db file gets created...
#V4.2:to produce a token per variant that incudes tissue_type(universal):pmids...
d2 = {
'condition1':{'all':'haematopoietic_and_lymphoid_tissue>HS1=acute_lymphoblastic_*','aml':'haematopoietic_and_lymphoid_tissue>HS1=acute_myeloid_*',
'mds':'haematopoietic and lymphoid tissue>HS1=myelodysplastic*',
'nhl':['haematopoietic_and_lymphoid_tissue>HS1=anaplastic large cell lymphoma',
'haematopoietic_and_lymphoid_tissue>HS1=angioimmunoblastic t cell lymphoma',
'haematopoietic_and_lymphoid_tissue>HS1=b cell prolymphocytic leukaemia',
'haematopoietic_and_lymphoid_tissue>HS1=burkitt lymphoma',
'haematopoietic_and_lymphoid_tissue>HS1=chronic lymphocytic leukaemia-small lymphocytic lymphoma',
'haematopoietic_and_lymphoid_tissue>HS1=diffuse large b cell lymphoma',
'haematopoietic_and_lymphoid_tissue>HS1=follicular lymphoma',
'haematopoietic_and_lymphoid_tissue>HS1=gamma-delta t cell lymphoma',
'haematopoietic_and_lymphoid_tissue>HS1=hairy cell leukaemia',
'haematopoietic_and_lymphoid_tissue>HS1=hepatosplenic t cell lymphoma',
'haematopoietic_and_lymphoid_tissue>HS1=lymphoplasmacytic lymphoma',
'haematopoietic_and_lymphoid_tissue>HS1=malt lymphoma,mantle cell lymphoma',
'haematopoietic_and_lymphoid_tissue>HS1=marginal zone lymphoma',
'haematopoietic_and_lymphoid_tissue>HS1=monoclonal gammopathy of undetermined significance',
'haematopoietic_and_lymphoid_tissue>HS1=nk-t cell lymphoma',
'haematopoietic_and_lymphoid_tissue>HS1=nk cell leukaemia',
'haematopoietic_and_lymphoid_tissue>HS1=peripheral t cell lymphoma unspecified',
'haematopoietic_and_lymphoid_tissue>HS1=plasma cell myeloma',
'haematopoietic_and_lymphoid_tissue>HS1=plasmacytoma',
'haematopoietic_and_lymphoid_tissue>HS1=primary central nervous system lymphoma',
'haematopoietic_and_lymphoid_tissue>HS1=primary effusion lymphoma',
'haematopoietic_and_lymphoid_tissue>HS1=t cell large granular lymphocytic leukaemia'],
'nsclc':'lung>HS1=non_small_cell_carcinoma','sclc':'lung>HS1=small_cell_carcinoma'},
'condition2':{'anal':'large_intestine>SS1=anus','rectal':'large_intestine>SS1=rectum','colon':'large_intestine>SS1=colon','bladder':'urinary_tract>SS1=bladder'},
'condition3':{'uterine':['cervix','endometrium'],'bone':'bone','breast':'breast','cns':'central_nervous_system','esophageal':'oesophagus','gastric':'stomach',
'ovarian':'ovary','cervical':'cervix','pancreas':'pancreas'},
'condition4':{'kaposi':'soft_tissue>PH=Kaposi_sarcoma','melanoma':'PH=malignant_melanoma','hnscc':'skin>SS1=head_neck>PH=carcinoma>HS1=Squamous cell carcinoma'}
}
#
try:
	if len(sys.argv) < 2:
		print ('USAGE:python CosmicMutantExportCensus.DDMMYYYY.tsv Tissue_Types.*.csv')
		exit(0)
except:
	pass
fH2 = open (sys.argv[2],'r')
lines2 = fH2.readlines()
tt={}
for line in lines2:
	#line  = line.decode('utf-7','ignore')
	line=line.replace('\n','').replace('+AD4-','>').replace('+ADs-',';').replace('+AC0-','-').replace('+ACo-','*').replace('+AD0-','=').replace('+AF8-','_').replace('+AF8AKg-','_*').replace('+ACY-','&')
	tokens =  line.split('\t')
	if tokens[0] == 'origin':
		continue
	else:
		tt[tokens[1]] = tokens[2]
#
fH = open (sys.argv[1],'r')
lines = fH.readlines()
D={}
for line in lines:
	#line  = line.decode('utf-7','ignore')
	tokens = line.split('\t')
	if tokens[0] =='Gene name':
		header = line
	else:
		var = tokens[0] +':'+tokens[1]+':'+tokens[17]
		tumor_site = tokens[7]+';'+tokens[8]+';'+tokens[9]+';'+tokens[10]+';'+tokens[11]+';'+tokens[12]+';'+tokens[13]+';'+tokens[14]#V4.2#

		D[var]={'Genename':'','HGNCid':'','UniversalTissueType':'NA','Primarysite':'','Sitesubtype1':'','Sitesubtype2':'','Sitesubtype3':'','Primaryhistology':'','Histologysubtype1':'','Histologysubtype2':'','Histologysubtype3':'','MutationID':'','AA':'','Mutationdescription':'','MutationGenomePosition':'','MutationSomaticStatus':'','TumourOrigin':'','PMID':'','Tier':'','token':{}}
	
for line in lines:
	#line  = line.decode('utf-7','ignore')
	tokens = line.split('\t')
	if tokens[0] =='Gene name':
		header = line
	else:
		var = tokens[0] +':'+tokens[1]+':'+tokens[17]
		#tumor_site = tokens[7]+';'+tokens[8]+';'+tokens[9]+';'+tokens[10]+';'+tokens[11]+';'+tokens[12]+';'+tokens[13]+';'+tokens[14]#V4.2#
		if var in D:
			D[var]['Genename'] += ';'+tokens[0];D[var]['HGNCid'] += ';'+tokens[3];D[var]['Primarysite'] += ';'+tokens[7];D[var]['Sitesubtype1'] += ';'+tokens[8];D[var]['Sitesubtype2'] += ';'+tokens[9];D[var]['Sitesubtype3'] += ';'+tokens[10];D[var]['Primaryhistology'] += ';'+tokens[11];D[var]['Histologysubtype1'] += ';'+tokens[12];D[var]['Histologysubtype2'] += ';'+tokens[13];D[var]['Histologysubtype3'] += ';'+tokens[14];D[var]['MutationID'] += ';'+tokens[16];D[var]['AA'] += ';'+tokens[18];D[var]['Mutationdescription'] += ';'+tokens[19];D[var]['MutationGenomePosition'] += ';'+tokens[23];D[var]['MutationSomaticStatus'] += ';'+tokens[29];D[var]['TumourOrigin'] += ';'+tokens[33];D[var]['PMID'] += ';'+tokens[30];D[var]['Tier'] += ';'+tokens[35].replace('\n','')
			#V4.2#
			#D[var]['token'].append(typ+':'+tokens[7]+';'+tokens[8]+';'+tokens[9]+';'+tokens[10]+';'+tokens[11]+';'+tokens[12]+';'+tokens[13]+';'+tokens[14]+':'+tokens[30])
			#V4.2#
			cond1 = tokens[7]+'>HS1='+tokens[12]
			cond2 = tokens[7]+'>SS1='+tokens[8] 
			cond3 = tokens[7]
			cond4_1 = tokens[7] +'>PH='+tokens[11]
			cond4_2 = tokens[11]
			cond4_3 = tokens[7] +'>SS1='+tokens[8]+'>PH='+tokens[11]+'>HS1='+tokens[12]
			for typ in d2['condition1']:
				if typ =='nhl':
					for t in d2['condition1'][typ]:
						if cond1 == t:
							D[var]['UniversalTissueType'] += ';' +typ
							#tumor_site += ':'+typ;D[var]['token'][tumor_site]=[]
							tumor_site = typ;D[var]['token'][tumor_site]=[]
				elif d2['condition1'][typ].replace('*','') in cond1:
					D[var]['UniversalTissueType'] += ';' +typ
					#tumor_site += ':'+typ;D[var]['token'][tumor_site]=[]
					tumor_site = typ;D[var]['token'][tumor_site]=[]
			for typ in d2['condition2']:
				if cond2 ==  d2['condition2'][typ]:
					D[var]['UniversalTissueType'] += ';' +typ
					#tumor_site += ':'+typ;D[var]['token'][tumor_site]=[]
					tumor_site = typ;D[var]['token'][tumor_site]=[]
			for typ in d2['condition3']:
				if typ == 'uterine':
					for t in d2['condition3'][typ]:
						if cond3 == t:
							D[var]['UniversalTissueType'] += ';' +typ
							#tumor_site += ':'+typ;D[var]['token'][tumor_site]=[]
							tumor_site = typ;D[var]['token'][tumor_site]=[]
				elif cond3 == d2['condition3'][typ]:
					D[var]['UniversalTissueType'] += ';' +typ
					#tumor_site += ':'+typ;D[var]['token'][tumor_site]=[]
					tumor_site = typ;D[var]['token'][tumor_site]=[]
			for typ in d2['condition4']:
				if d2['condition4'][typ] == cond4_1:
					D[var]['UniversalTissueType'] += ';' +typ
					#tumor_site += ':'+typ;D[var]['token'][tumor_site]=[]
					tumor_site = typ;D[var]['token'][tumor_site]=[]
				elif d2['condition4'][typ] == cond4_2:
					D[var]['UniversalTissueType'] += ';' +typ
					#tumor_site += ':'+typ;D[var]['token'][tumor_site]=[]
					tumor_site = typ;D[var]['token'][tumor_site]=[]
				elif d2['condition4'][typ] == cond4_3:
					D[var]['UniversalTissueType'] += ';' +typ
					#tumor_site += ':'+typ;D[var]['token'][tumor_site]=[]
					tumor_site = typ;D[var]['token'][tumor_site]=[]
		else:
			print ('WARNING:Could not find this variant:{}'.format(var))
			exit(0)
#4.2
for line in lines:
	#line  = line.decode('utf-7','ignore')
	tokens = line.split('\t')
	if tokens[0] =='Gene name':
		header = line
	else:
		var = tokens[0] +':'+tokens[1]+':'+tokens[17]
		#tumor_site = tokens[7]+';'+tokens[8]+';'+tokens[9]+';'+tokens[10]+';'+tokens[11]+';'+tokens[12]+';'+tokens[13]+';'+tokens[14]#V4.2#
		if var in D:
			cond1 = tokens[7]+'>HS1='+tokens[12]
			cond2 = tokens[7]+'>SS1='+tokens[8] 
			cond3 = tokens[7]
			cond4_1 = tokens[7] +'>PH='+tokens[11]
			cond4_2 = tokens[11]
			cond4_3 = tokens[7] +'>SS1='+tokens[8]+'>PH='+tokens[11]+'>HS1='+tokens[12]
			for typ in d2['condition1']:
				if typ =='nhl':
					for t in d2['condition1'][typ]:
						if cond1 == t:
							#tumor_site += ':'+typ;D[var]['token'][tumor_site].append(tokens[30])
							tumor_site = typ;D[var]['token'][tumor_site].append(tokens[30])
				elif d2['condition1'][typ].replace('*','') in cond1:
					#tumor_site += ':'+typ;D[var]['token'][tumor_site].append(tokens[30])
					tumor_site = typ;D[var]['token'][tumor_site].append(tokens[30])
			for typ in d2['condition2']:
				if cond2 ==  d2['condition2'][typ]:
					#tumor_site += ':'+typ;D[var]['token'][tumor_site].append(tokens[30])
					tumor_site = typ;D[var]['token'][tumor_site].append(tokens[30])
			for typ in d2['condition3']:
				if typ == 'uterine':
					for t in d2['condition3'][typ]:
						if cond3 == t:
							#tumor_site += ':'+typ;D[var]['token'][tumor_site].append(tokens[30])
							tumor_site = typ;D[var]['token'][tumor_site].append(tokens[30])
				elif cond3 == d2['condition3'][typ]:
					#tumor_site += ':'+typ;D[var]['token'][tumor_site].append(tokens[30])
					tumor_site = typ;D[var]['token'][tumor_site].append(tokens[30])
			for typ in d2['condition4']:
				if d2['condition4'][typ] == cond4_1:
					#tumor_site += ':'+typ;D[var]['token'][tumor_site].append(tokens[30])
					tumor_site = typ;D[var]['token'][tumor_site].append(tokens[30])
				elif d2['condition4'][typ] == cond4_2:
					#tumor_site += ':'+typ;D[var]['token'][tumor_site].append(tokens[30])
					tumor_site = typ;D[var]['token'][tumor_site].append(tokens[30])
				elif d2['condition4'][typ] == cond4_3:
					#tumor_site += ':'+typ;D[var]['token'][tumor_site].append(tokens[30])
					tumor_site = typ;D[var]['token'][tumor_site].append(tokens[30])
		else:
			print ('WARNING:Could not find this variant:{}'.format(var))
			exit(0)
#4.2
print (len(D))
names = sys.argv[1].split('.')
ofile = 'COSMIC.db.'+names[-2]+'.v4.csv'
fOut = open (ofile,'w')
fOut.write('variant\tGene name\tHGNC ID\tUniversal type\tPrimary site\tSite subtype1\tSite subtype2\tSite subtype3\tPrimary histology\tHistology subtype1\tHistology subtype2\tHistology subtype3\tMutation ID\tMutation AA\tMutation description\tMutation genome position\tMutation somatic status\tTumour origin\tPMID\tTier\ttoken\n')
I=1
for var in D:
	fOut.write('{}\t'.format(var))
	for k in D[var]:
		if k == 'token':#'Tier':
			#tokens = D[var][k].split(';')
			#tokens = list(set(tokens))
			for ts in D[var][k]:
				pmids = list(set(D[var][k][ts]))
				s = ts +':'+ ','.join(pmids)
				fOut.write('{}#'.format(s))
			fOut.write('\n')	
		else:
			tokens = D[var][k].split(';')
			tokens = list(set(tokens))
			s = ';'.join(tokens)
			fOut.write('{}\t'.format(s))
	#print ('NOTICE:Done {} of {}'.format(I,len(D)));I+=1
print ('NOTICE:Done writing {} variants into {}...'.format(len(D),ofile))







