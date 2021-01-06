#!/usr/bin/env python3
import os
import sys
import collections
#This is to classify somatic variants given their name & tumor type (gene:ENST_id:cDNA ttype)...
#v3(SV-Class):To increase scripts dynamic...
#3.3:To correct the LOGIC of wrongly-classified variants that should be in class III not I because B level of evidence is tumor site specific based on the guidelines examples...
#3.4: Remove not necessary (unused) db...
try:
	if len(sys.argv) < 2:
		print('USAGE:python input_vairants_list_files.txt')
		exit(0)
except:
	pass
#functions...
#1.Reading dbs...#MCG already covered in PG!
def finalDB(path):
	os.chdir(path)
	files = os.listdir(path)
	#files = path.split(',')
	COSMIC={};cBioPortal={};OncoKB={};OncoKB_unspecified={}
	CENSUS={};types_universal={};types_original={}
	PG={}
	for f in files:
		#dirc = f.split('/')
		name = f.split('.')
		if name[0] =='COSMIC':
			fH = open (f,'r')
			lines = fH.readlines()
			del lines[0]
			for line in lines:
				line = line.replace('\n','')
				tokens = line.split('\t')
				COSMIC[tokens[0]] =[]
				if tokens[-1] != '':#to avoid variants with unclear site(there are many!)....
					parts = tokens[-1].split('#')
					for y in range(0,len(parts)-1,++1):
						info = parts[y].split(':')
						pmids = info[1].split(',')
						pmids = list(set(pmids))
						COSMIC[tokens[0]].append([info[0],len(pmids)])
				#pmids = tokens[-2].split(';')
				#pmid = 0 
				#for i in pmids:
				#	if i =='':
				#		continue
				#	else:
				#		pmid+=1
				#COSMIC[tokens[0]] = [tokens[3],pmid]
			fH.close()
		elif name[0] == 'cBioPortal':
			fH = open (f,'r')
			lines = fH.readlines()
			del lines[0]
			for line in lines:
				line = line.replace('\n','')
				tokens = line.split('\t')
				parts = tokens[0].split(':')
				try:
					chuncks = parts[1].split('.')
					var = parts[0]+':'+chuncks[0]+':'+parts[2]
					cBioPortal[var] =[]
					refs = tokens[-1].split('#')
					for r in range(0,len(refs)-1,++1):
						info = refs[r].split(':')
						pmids = info[1].split(';')
						pmids = list(set(pmids))
						cBioPortal[var].append([info[0],info[2],len(pmids)])#NEWCODE...
					#pmids = tokens[-2].split(';')
					#pmid = 0
					#for i in pmids:
					#	if i =='NA':
					#		continue
					#	else:
					#		pmid+=1
							   #group      type      pmids_count
					#cBioPortal[var] = [tokens[-3],tokens[-1],pmid]#tokens[-3] to measure volume/well-powered study...
				except:
					pass
			fH.close()
		
		elif name[0] == 'OncoKB':
			fH = open (f,'r')
			lines = fH.readlines()
			del lines[0]
			for line in lines:
				tokens = line.split('\t')
				OncoKB[tokens[0]]=[]
			for line in lines:
				tokens = line.split('\t')
				gene = tokens[0];del tokens[0]
				OncoKB[gene].append(tokens)
			fH.close()
		elif name[0] =='OncoKB_unspecified-alterations':
			fH = open (f,'r')
			lines = fH.readlines()
			del lines[0]
			for line in lines:
				line = line.replace('\n','')
				tokens  = line.split('\t')
				gene = tokens[0];del tokens[0]
				OncoKB_unspecified[gene]=tokens
			fH.close()
		
		elif name[0] =='CENSUS':
			fH = open (f,'r')
			lines = fH.readlines()
			del lines[0]
			for line in lines:
				line = line.replace('\n','')
				tokens = line.split('\t')
				gene = tokens[0];del tokens[0];del tokens[0]
				CENSUS[gene] = {'all':tokens[0], 'aml':tokens[1], 'anal':tokens[2], 'rectal':tokens[3], 'colon':tokens[4], 'bladder':tokens[5], 'uterine':tokens[6], 'bone':tokens[7], 'breast':tokens[8], 'cns':tokens[9], 'esophageal':tokens[10], 'gastric':tokens[11], 'kaposi':tokens[12], 'melanoma':tokens[13], 'ovarian':tokens[14], 'cervical':tokens[15], 'mds':tokens[16], 'nhl':tokens[17], 'pancreas':tokens[18], 'nsclc':tokens[19], 'sclc':tokens[20], 'hnscc':tokens[21]}
			fH.close()
		
		elif name[0] =='Tumor_Type':
			#two dictionaries will be returned...
			fH = open (f,'r')
			lines = fH.readlines()
			del lines[0]
			for line in lines:
				line = line.replace('\n','').replace('+AD4-','>').replace('+ADs-',';').replace('+AC0-','-').replace('+ACo-','*').replace('+AD0-','=').replace('+AF8-','_').replace('+AF8AKg-','_*').replace('+ACY-','&')
				tokens = line.split('\t')
				universal = tokens[1];del tokens[1]
				types_universal[universal]=tokens
				types_original[tokens[0].lower()] = []
			for line in lines:
				line = line.replace('\n','').replace('+AD4-','>').replace('+ADs-',';').replace('+AC0-','-').replace('+ACo-','*').replace('+AD0-','=').replace('+AF8-','_').replace('+AF8AKg-','_*').replace('+ACY-','&')
				tokens = line.split('\t')
				types_original[tokens[0].lower()].append(tokens[1])
			fH.close()
		
		elif name[0] == 'PGs_all':
			fH = open (f,'r')
			lines = fH.readlines()
			del lines[0]
			for line in lines:
				line = line.replace('+AD4-','>').replace('+ADs-',';').replace('+AC0-','-').replace('+ACo-','*').replace('+ACY-','&')
				tokens = line.split('\t')
				if tokens[0] =='driver_mutations_not_found_in_db':
					continue
				else:
					genes = tokens[1].split(';')
					for gene in genes:
						if gene != '':
							PG[gene]=[]
			for line in lines:
				line = line.replace('+AD4-','>').replace('+ADs-',';').replace('+AC0-','-').replace('+ACo-','*').replace('+ACY-','&')
				tokens = line.split('\t')
				if tokens[0] =='driver_mutations_not_found_in_db':
					continue
				else:
					genes = tokens[1].split(';')
					del tokens[1]
					for gene in genes:
						if gene != '':
							PG[gene].append(tokens)
			fH.close()
		else:
			continue
		

	return COSMIC,cBioPortal,OncoKB,OncoKB_unspecified,CENSUS,types_universal,types_original,PG
def Infile(f):
	fH = open (f,'r')
	lines = fH.readlines()
	q={}
	for line in lines:
		tokens = line.split('\t')
		parts = tokens[0].split(':')
		gene = parts[0]
		q[gene]=[]
	for line in lines:
		tokens = line.split('\t')
		parts = tokens[0].split(':')
		gene = parts[0]
		typ = tokens[2].replace('\n','')
		del tokens[2]
		q[gene].append(tokens)
		
	fH.close()
	return q,typ
path = sys.argv[2].replace('\n','')#'/path/finalDB/'

COSMIC,cBioPortal,OncoKB,OncoKB_unspecified,CENSUS,types_universal,types_original,PG = finalDB(path)

print('DONE:DB READING...')
ipath = sys.argv[1]#'/path/Input/tumor_site.query_sample.txt' 
try:
	genes, typ = Infile(ipath)
	print('DONE:INFILE READING...')
except Exception as e:
	print(e)
        
#avoid to get over variants in OncoKB that either not captured by DNA-Seq, need require a VCF file or not found in intOgen... 
avoid = ['Amplification','BCR-ABL1 Fusion','COL1A1-PDGFB Fusion','EWSR1-FLI1 Fusion','Exon 14 Deletion','FIP1L1-PDGFRA Fusion','Fusions','Internal tandem duplication','Kinase Domain Duplication','Microsatellite Instability-High','PCM1-JAK2 Fusion','Wildtype']
output={};done=[]
for gene in genes:
	for var in genes[gene]:
		tokens = var[0].split(':')
		output[var[0]] = {'gene':gene,'ens':tokens[1],'cDNA':tokens[2],'aa':var[1],'submitted_type':typ,'tier':'','level':'','bio_type':set(),'NCCN_drug':'','NCCN_response':'','NCCN_condition':'','NCCN_panel':'','FDA_for_other':set(),'FDA_for_same':set(),'investigational_for_same':set(),'preclinical_generic':set(),'preclinical_for_same':set(),'ct_title':set(),'status':set(),'location':set(),'e_source':set()}
L=[];Qu=[]
track=0
for gene in genes:
	#level A...
	if gene in PG:
		for var in genes[gene]:
			for i in PG[gene]:
				if typ in i[-3]:#get tissue type...
					var[1] = var[1].replace('p.','')
					if i[0][-1] == '*':#remove only at the end...
						i[0] = i[0].replace('*','')
					tmp = var[1].split('_')#to remove del/ins...
					aa ='na'
					if len(tmp) == 1:
						aa = tmp[0]
					if var[1] in i[0]:#search by aa...
						done.append(var[0])
						#output[var[0]]['tier'] = 'I'
						output[var[0]]['tier'] = 'I';output[var[0]]['level']='A';output[var[0]]['bio_type'].add(i[1])
						output[var[0]]['NCCN_drug'] =  i[2]
						output[var[0]]['NCCN_response'] =  i[3]
						output[var[0]]['NCCN_condition'] =  i[5]
						output[var[0]]['NCCN_panel'] =   i[6]
						output[var[0]]['e_source'].add('NCCN')
					elif i[0] in aa:#search by modified aa (cover both IARC & BRCAEx)...
						done.append(var[0])
						#output[var[0]]['tier'] ='I'
						output[var[0]]['tier'] = 'I';output[var[0]]['level']='A';output[var[0]]['bio_type'].add(i[1])
						output[var[0]]['NCCN_drug'] =  i[2]
						output[var[0]]['NCCN_response'] =  i[3]
						output[var[0]]['NCCN_condition'] =  i[5]
						output[var[0]]['NCCN_panel'] =   i[6]
						output[var[0]]['e_source'].add('NCCN')
					elif var[0] in i[0]:#search by entire variant...
						done.append(var[0])
						#output[var[0]]['tier']='I'
						output[var[0]]['tier'] = 'I';output[var[0]]['level']='A';output[var[0]]['bio_type'].add(i[1])
						output[var[0]]['NCCN_drug'] =  i[2]
						output[var[0]]['NCCN_response'] =  i[3]
						output[var[0]]['NCCN_condition'] =  i[5]
						output[var[0]]['NCCN_panel'] =   i[6]
						output[var[0]]['e_source'].add('NCCN')
	if gene in OncoKB:
		for i in OncoKB[gene]:
			if i[-3] in types_universal[typ][-4]:#get tissue type...
				if i[-2] == '1':#get FDA approved only...
					for var in genes[gene]:
						var[1] = var[1].replace('p.','')
						if var[1] == i[1]:
							done.append(var[0])#
							output[var[0]]['tier'] = 'I';output[var[0]]['level']='A';output[var[0]]['bio_type'].add('therapeutic')
							output[var[0]]['FDA_for_same'].add(i[-1])
							output[var[0]]['e_source'].add('OncoKB')
						elif gene in OncoKB_unspecified:
							if i[1] not in avoid:#11.2.2020						
								if var[0] in '\t'.join(OncoKB_unspecified[gene]):#get un-specified variants...
									done.append(var[0])
									output[var[0]]['tier'] = 'I';output[var[0]]['level']='A';output[var[0]]['bio_type'].add('therapeutic')
									output[var[0]]['FDA_for_same'].add(i[-1])
									output[var[0]]['e_source'].add('OncoKB')	
	#level B & C...
	for var in genes[gene]:
		Qu.append(var[0])
		if var[0] in done:
			continue
		else:
			try:
				for varInfo in cBioPortal[var[0]]:
					if 'PANCAN' in varInfo[0]:#well-powered studies...
						if gene in CENSUS:
							if int(CENSUS[gene][typ]) == 1:#with COSMIC census (tier I or II)...#11.2.2020
								done.append(var[0])
								output[var[0]]['tier'] = 'I';output[var[0]]['level']='B';output[var[0]]['bio_type'].add('diagnostic/prognostic')
								output[var[0]]['e_source'].add('PAN_with_COSMIC_Consensus')#11.2.2020
								#get FDA approved drugs and investigational ones JUST FOR INFORMATION WITHOUT AFFECTING CLASSIFICATION...
								if gene in OncoKB:
									for i in OncoKB[gene]:
										if i[-3] in types_universal[typ][-4]:#get tissue type...
											if i[-2] == '3A' or i[-2] == '3B':#get investigational drug only...
												for var in genes[gene]:
													var[1] = var[1].replace('p.','')
													if var[1] == i[1]:
														done.append(var[0])
														output[var[0]]['investigational_for_same'].add(i[-1])
													elif var[0] in '\t'.join(OncoKB_unspecified[gene]):#get un-specified variants...
														done.append(var[0])
														output[var[0]]['investigational_for_same'].add(i[-1])
										else:#get FDA approved only regardless of tissue type...
											if i[-2] == '1':
												for var in genes[gene]:
													var[1] = var[1].replace('p.','')
													if var[1] == i[1]:
														done.append(var[0])
														output[var[0]]['FDA_for_other'].add(i[-1])
													elif gene in OncoKB_unspecified:
														if i[1] not in avoid:
															if var[0] in '\t'.join(OncoKB_unspecified[gene]):#get un-specified variants...
																done.append(var[0])
																output[var[0]]['FDA_for_other'].add(i[-1])
								##########
							elif int(CENSUS[gene][typ]) == 2:#with CCGD census...#11.2.2020_all_statement_way_down...
								done.append(var[0])
								output[var[0]]['tier'] = 'II';output[var[0]]['level']='C';output[var[0]]['bio_type'].add('diagnostic/prognostic')
								output[var[0]]['e_source'].add('PAN_with_CCGD_Consensus')#1.2.2020
								#get FDA approved drugs and investigational ones JUST FOR INFORMATION WITHOUT AFFECTING CLASSIFICATION...
								if gene in OncoKB:
									for i in OncoKB[gene]:
										if i[-3] in types_universal[typ][-4]:#get tissue type...
											if i[-2] == '3A' or i[-2] == '3B':#get investigational drug only...
												for var in genes[gene]:
													var[1] = var[1].replace('p.','')
													if var[1] == i[1]:
														done.append(var[0])
														output[var[0]]['investigational_for_same'].add(i[-1])
													elif var[0] in '\t'.join(OncoKB_unspecified[gene]):#get un-specified variants...
														done.append(var[0])
														output[var[0]]['investigational_for_same'].add(i[-1])
										else:#get FDA approved only regardless of tissue type...
											if i[-2] == '1':
												for var in genes[gene]:
													var[1] = var[1].replace('p.','')
													if var[1] == i[1]:
														done.append(var[0])
														output[var[0]]['FDA_for_other'].add(i[-1])
													elif gene in OncoKB_unspecified:
														if i[1] not in avoid:
															if var[0] in '\t'.join(OncoKB_unspecified[gene]):#get un-specified variants...
																done.append(var[0])
																output[var[0]]['FDA_for_other'].add(i[-1])
			
					elif 'PANCAN' not in varInfo[0]:
						if varInfo[2] > 5:#another census layer#11.2.2020(moved)#multiple small studies...



							if gene in CENSUS:
								if int(CENSUS[gene][typ]) == 1 or int(CENSUS[gene][typ]) == 2:#with any census#11.2.2020
										done.append(var[0])		
										output[var[0]]['tier'] = 'II';output[var[0]]['level']='C';output[var[0]]['bio_type'].add('diagnostic/prognostic')					
										output[var[0]]['e_source'].add('Multiple_with_some_consensus')
						elif varInfo[2] < 5:#case reports...
							
							if typ == varInfo[1]:#to get related tissue type...
								done.append(var[0])
								output[var[0]]['tier'] = 'II';output[var[0]]['level']='D';output[var[0]]['bio_type'].add('diagnostic/prognostic')									
								output[var[0]]['e_source'].add('cBioPortal')
				#output[var[0]]['e_source'].add('cBioPortal')
							
			except:#Class II,level C...  				
				if gene in OncoKB:
					for i in OncoKB[gene]:
						if i[-3] in types_universal[typ][-4]:#get tissue type...
							if i[-2] == '3A' or i[-2] == '3B':#get investigational drug only...
								for var in genes[gene]:
									var[1] = var[1].replace('p.','')
									if var[1] == i[1]:
										done.append(var[0])
										if output[var[0]]['tier'] != 'I':
											output[var[0]]['tier'] = 'II';output[var[0]]['level']='C';output[var[0]]['bio_type'].add('therapeutic')
											output[var[0]]['investigational_for_same'].add(i[-1])
											output[var[0]]['e_source'].add('OncoKB')
										elif output[var[0]]['tier'] == 'I':
											output[var[0]]['investigational_for_same'].add(i[-1])
									elif var[0] in '\t'.join(OncoKB_unspecified[gene]):#get un-specified variants...
										done.append(var[0])
										if output[var[0]]['tier'] != 'I':
											output[var[0]]['tier'] = 'II';output[var[0]]['level']='C';output[var[0]]['bio_type'].add('therapeutic')
											output[var[0]]['investigational_for_same'].add(i[-1])
											output[var[0]]['e_source'].add('OncoKB')
										elif output[var[0]]['tier'] =='I':
											output[var[0]]['investigational_for_same'].add(i[-1])
						else:#get FDA approved only regardless of tissue type...
							if i[-2] == '1':
								for var in genes[gene]:
									var[1] = var[1].replace('p.','')
									if var[1] == i[1]:
										done.append(var[0])
										if output[var[0]]['tier'] != 'I':
											output[var[0]]['tier'] = 'II';output[var[0]]['level']='C';output[var[0]]['bio_type'].add('therapeutic')
											output[var[0]]['FDA_for_other'].add(i[-1]+':'+i[-3])
											output[var[0]]['e_source'].add('OncoKB')
										elif output[var[0]]['tier'] =='I':
											output[var[0]]['FDA_for_other'].add(i[-1]+':'+i[-3])
									elif gene in OncoKB_unspecified:
										if i[1] not in avoid:
											if var[0] in '\t'.join(OncoKB_unspecified[gene]):#get un-specified variants...
												done.append(var[0])
												if output[var[0]]['tier'] != 'I':
													output[var[0]]['tier'] = 'II';output[var[0]]['level']='C';output[var[0]]['bio_type'].add('therapeutic')
													output[var[0]]['FDA_for_other'].add(i[-1]+':'+i[-3])
													output[var[0]]['e_source'].add('OncoKB')
												elif output[var[0]]['tier'] == 'I':
													output[var[0]]['FDA_for_other'].add(i[-1]+':'+i[-3])
				
				if var[0] in COSMIC:
					for varInfo in COSMIC[var[0]]:
						if varInfo[1] > 5:#mutiple small studies...
							if gene in CENSUS:
								if int(CENSUS[gene][typ]) == 1 or int(CENSUS[gene][typ]) == 2:#with any census#11.2.2020
									
									if typ == varInfo[0]:#to get tissue type specific...	
										done.append(var[0])
										if output[var[0]]['tier'] != 'I':
											output[var[0]]['tier'] = 'II';output[var[0]]['level']='C';output[var[0]]['bio_type'].add('diagnostic/prognostic')												
											output[var[0]]['e_source'].add('COSMIC')
										elif output[var[0]]['tier'] == 'II' and output[var[0]]['level'] == 'C':
											output[var[0]]['bio_type'].add('diagnostic/prognostic')
											output[var[0]]['e_source'].add('COSMIC')
						elif COSMIC[var[0]][1] < 5:#Class D assuming census score of 0 or 1...
							
							if typ == varInfo[0]:#to get tissue type specific...	
								done.append(var[0])
								if output[var[0]]['tier'] != 'I':
									output[var[0]]['tier'] = 'II';output[var[0]]['level']='D';output[var[0]]['bio_type'].add('diagnostic/prognostic')										
									output[var[0]]['e_source'].add('COSMIC')
								elif output[var[0]]['tier'] == 'II' and output[var[0]]['level'] == 'D':
									output[var[0]]['bio_type'].add('diagnostic/prognostic')
									output[var[0]]['e_source'].add('COSMIC')
						#Will skip intOgen because there is no number of studies and it has not been updated in a long time....
		
	#level D....
	for var in genes[gene]:
		if var[0] in done:
			continue
		else:
			if gene in OncoKB:
				for i in OncoKB[gene]:
					if i[-3] in types_universal[typ][-4] :#get tissue type ...
						if i[-2] == '4':
							for var in genes[gene]:
								var[1] = var[1].replace('p.','')
								if var[1] == i[1]:
									done.append(var[0])
									if output[var[0]]['tier'] =='':
										output[var[0]]['tier'] = 'II';output[var[0]]['level']='D';output[var[0]]['bio_type'].add('therapeutic')
										output[var[0]]['preclinical_for_same'].add(i[-1])
										output[var[0]]['e_source'].add('OncoKB')
									elif output[var[0]]['tier'] !='':
										output[var[0]]['preclinical_for_same'].add(i[-1])
								elif gene in OncoKB_unspecified:
									if i[1] not in avoid:
										if var[0] in '\t'.join(OncoKB_unspecified[gene]):#get un-specified variants...
											done.append(var[0])
											if output[var[0]]['tier'] =='':
												output[var[0]]['tier'] = 'II';output[var[0]]['level']='D';output[var[0]]['bio_type'].add('therapeutic')
												output[var[0]]['preclinical_for_same'].add(i[-1])
												output[var[0]]['e_source'].add('OncoKB')
											elif output[var[0]]['tier'] !='':
												output[var[0]]['preclinical_for_same'].add(i[-1])
					elif i[-3] == 'All Solid Tumors' or i[-3] == 'All Tumors':#generic treatments...
						if i[-2] == '4':
							for var in genes[gene]:
								var[1] = var[1].replace('p.','')
								if var[1] == i[1]:
									done.append(var[0])
									if output[var[0]]['tier'] =='':
										output[var[0]]['tier'] = 'II';output[var[0]]['level']='D';output[var[0]]['bio_type'].add('therapeutic')
										output[var[0]]['preclinical_generic'].add(i[-1]+':'+i[-3])
										output[var[0]]['e_source'].add('OncoKB')
									elif output[var[0]]['tier'] !='':
										output[var[0]]['preclinical_generic'].add(i[-1]+':'+i[-3])
								elif gene in OncoKB_unspecified:
									if i[1] not in avoid:
										if var[0] in '\t'.join(OncoKB_unspecified[gene]):#get un-specified variants...
											done.append(var[0])
											if output[var[0]]['tier'] =='':
												output[var[0]]['tier'] = 'II';output[var[0]]['level']='D';output[var[0]]['bio_type'].add('therapeutic')
												output[var[0]]['preclinical_generic'].add(i[-1]+':'+i[-3])
												output[var[0]]['e_source'].add('OncoKB')
											elif output[var[0]]['tier'] !='':
												output[var[0]]['preclinical_generic'].add(i[-1]+':'+i[-3])
						
	#Class 3
	for var in genes[gene]:
		if var[0] in done:
				continue
		else:
			output[var[0]]['tier'] = 'III'#This is not "classifed" and is ultra rare or absent from pupblic database...
			
					
	for var in genes[gene]:
		if var[0] in done:
			continue
		else:
			L.append(var[0])

L= list(set(L))#L == left...
done = list(set(done))
Qu = list(set(Qu))#Qu == Quality...	
#
print(len(done));print(len(L))
print(len(done) + len(L))
print(len(Qu))
###Checking on the output for selected variants to cover all conditions...
#print(output['KIT:ENST00000288135:c.1924A>G'])#I/A
#print(output['BRAF:ENST00000288602:c.1799T>A'])#I/A
#print(output['GATAD2B:ENST00000368655:c.124C>T'])#I/B
#print(output['BRF2:ENST00000220659:c.1256C>T'])#I/B
#print(output['CDC42EP1:ENST00000249014:c.762_782delGCCTGCTGCAAACCCCTCAGC'])#II/C
#print(output['CRYM:ENST00000219599:c.749C>T'])#II/D
#print(output['ERBB2:ENST00000269571:c.929C>T'])#II/C
#print(output['EPCAM:ENST00000263735:c.475A>T'])#III/E
#print(output['MTOR:ENST00000361445:c.6514T>A'])#II/D generic treatment

##formatting & writing the output...

ofile = sys.argv[1].replace('.txt','.classified.txt')
#ofile = ('/project/k1376/tmp.classified.txt')
print('..{}..'.format(ofile))
fOut = open (ofile,'w')
fOut.write('gene\tens\tcDNA\taa\ttier\tlevel\te_source\tbio_type\tsubmitted_type\tNCCN_drug\tNCCN_response\tNCCN_condition\tNCCN_panel\tFDA_for_same\tFDA_for_other\tinvestigational_for_same\tct_title\tstatus\tlocation\tpreclinical_for_same\tpreclinical_generic\n')
for var in output:
	for x in range(21):#to get information sorted...
		if x == 0:
			fOut.write('{}\t'.format(output[var]['gene']))
		elif x == 1:
			fOut.write('{}\t'.format(output[var]['ens']))
		elif x == 2:	
			fOut.write('{}\t'.format(output[var]['cDNA']))
		elif x == 3:	
			fOut.write('{}\t'.format(output[var]['aa']))
		elif x == 4:	
			fOut.write('{}\t'.format(output[var]['tier']))
		elif x == 5:	
			fOut.write('{}\t'.format(output[var]['level']))
		elif x == 6:
			out = ';'.join(output[var]['e_source'])	
			fOut.write('{}\t'.format(out))
		elif x == 7:
			out = ';'.join(output[var]['bio_type'])
			fOut.write('{}\t'.format(out))
		elif x == 8:	
			fOut.write('{}\t'.format(output[var]['submitted_type']))
		elif x == 9:	
			fOut.write('{}\t'.format(output[var]['NCCN_drug']))
		elif x == 10:	
			fOut.write('{}\t'.format(output[var]['NCCN_response']))
		elif x == 11:	
			fOut.write('{}\t'.format(output[var]['NCCN_condition']))
		elif x == 12:	
			fOut.write('{}\t'.format(output[var]['NCCN_panel'].replace('\n','')))
		elif x == 13:
			out = ';'.join(output[var]['FDA_for_same']).replace('\n','')	
			fOut.write('{}\t'.format(out))
		elif x == 14:
			out = ';'.join(output[var]['FDA_for_other']).replace('\n','')	
			fOut.write('{}\t'.format(out))
		elif x == 15:
			out = ';'.join(output[var]['investigational_for_same']).replace('\n','')	
			fOut.write('{}\t'.format(out))
		elif x == 16:
			out = ';'.join(output[var]['ct_title'])	
			fOut.write('{}\t'.format(out))
		elif x == 17:
			out = ';'.join(output[var]['status'])	
			fOut.write('{}\t'.format(out))
		elif x == 18:
			out = ':location:'.join(output[var]['location'])	
			fOut.write('{}\t'.format(out))
		elif x == 19:
			out = ';'.join(output[var]['preclinical_for_same']).replace('\n','')	
			fOut.write('{}\t'.format(out))
		elif x == 20:
			out = ';'.join(output[var]['preclinical_generic']).replace('\n','')	
			fOut.write('{}\n'.format(out))
		
