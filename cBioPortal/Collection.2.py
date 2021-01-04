#!/usr/bin/env python3
import os
import sys
#This is to map clinical samples & mutations extended files into dictionary of directories names...
try:
	if len(sys.argv) < 2:
		print ('USAGE:python TumorTypes.V2.csv')
		exit(0)
except:
	pass
fH = open (sys.argv[1],'r')
lines = fH.readlines()
D={}#directories dictionary...
cBio={'MExtended':{},'CSamples':{},'MStudy':{}}
del lines[0]
for line in lines:
	tokens = line.split('\t')
	if tokens[0] !='':
		D[tokens[0].replace('+AF8-','_')]=tokens[5].replace('\n','')
		cBio['MExtended'][tokens[0].replace('+AF8-','_')]=[]
print ('NOTICE:Found {} directories in {}...'.format(len(D),sys.argv[1]))
#looping over directories...
path ='/media/turki/HD-SL2/SV-db/cBioPortal_all/datahub/public/'
Dirs = os.listdir(path)
headers=['Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification','Variant_Type','Reference_Allele','Tumor_Seq_Allele2','HGVSc','HGVSp_Short']#,'EXON','Codons','BIOTYPE','DOMAINS']
print ('>>>{}<<<'.format(len(headers)))
ProbFiles=['nhl_bcgsc_2013','brca_broad','skcm_broad_dfarber','brca_mbcproject_wagle_2017','breast_msk_2018','brca_tcga','paad_qcmg_uq_2016','paad_icgc','hnc_mskcc_2016','skcm_broad_brafresist_2012','skcm_yale','hnsc_tcga_pub','brca_bccrc','es_dfarber_broad_2014']
info = ['name','description','pmid','groups']
for Dir in D:
	if Dir in Dirs:
		try:#mutations extended file...
			fH1 = open ('{}{}/data_mutations_extended.txt'.format(path,Dir),'r')
			lines1 = fH1.readlines()
			sampleID='NA';tmpD={}
			if Dir == 'brca_sanger':
				l=[]
				for line in lines1:
					tokens = line.split('\t')
					if tokens[0] == 'Hugo_Symbol':
						sampleID=tokens[15]
						for x in range(0,len(tokens),++1):
							if tokens[x].replace('\n','') in headers:
								#l.append(tokens[x].replace('\n',''))
								l.append('{}:{}'.format(tokens[x].replace('\n',''),x))
					else:
						l.sort()
						if len(l) ==len(headers):
							tmpD[tokens[15]]=[]
							tmpL=[]
							for i in l:
								chunks = i.split(':')
								tmpL.append(tokens[int(chunks[1])].replace('\n',''))
						else:#aviod problematic lines...
							continue
				l=[]
				for line in lines1:
					tokens = line.split('\t')
					if tokens[0] == 'Hugo_Symbol':
						sampleID=tokens[15]
						for x in range(0,len(tokens),++1):
							if tokens[x].replace('\n','') in headers:
								#l.append(tokens[x].replace('\n',''))
								l.append('{}:{}'.format(tokens[x].replace('\n',''),x))
					else:
						l.sort()
						if len(l) ==len(headers):
							tmpL=[]
							for i in l:
								chunks = i.split(':')
								tmpL.append(tokens[int(chunks[1])].replace('\n',''))
							tmp = ','.join(tmpL)
							tmpD[tokens[15]].append(tmp)
						else:#aviod problematic lines...
							continue
#				print ('{}:{}'.format(Dir,sampleID))
			else:
				l=[]
				for line in lines1:
					tokens = line.split('\t')
					if tokens[0] == 'Hugo_Symbol':
						sampleID=tokens[16]
						for x in range(0,len(tokens),++1):
							if tokens[x] in headers:
								l.append('{}:{}'.format(tokens[x],x))
					else:
						l.sort()
						if len(l) ==len(headers):
							tmpD[tokens[16]]=[]
							tmpL=[]
							for i in l:
								chunks = i.split(':')
								tmpL.append(tokens[int(chunks[1])].replace('\n',''))
						else:#aviod problematic lines...
							continue
				l=[]
				for line in lines1:
					tokens = line.split('\t')
					if tokens[0] == 'Hugo_Symbol':
						sampleID=tokens[16]
						for x in range(0,len(tokens),++1):
							if tokens[x] in headers:
								l.append('{}:{}'.format(tokens[x],x))
					else:
						l.sort()
						if len(l) ==len(headers):
							tmpL=[]
							for i in l:
								chunks = i.split(':')
								tmpL.append(tokens[int(chunks[1])].replace('\n',''))
							tmp = ','.join(tmpL)
							tmpD[tokens[16]].append(tmp)
						else:#aviod problematic lines...
							continue
#				print ('{}:{}'.format(Dir,sampleID))
			cBio['MExtended'][Dir]=tmpD
			
		except: 
			print ('Error:Could not find data_mutations_extended.txt in {}'.format(Dir))
			exit(0)
		try:#clinical samples files...
			if Dir in ProbFiles:#get filtered files...
				fH2 = open ('{}{}/data_clinical_samples.filtered.txt'.format(path,Dir),'r')#;print ('######{}'.format(Dir))
			else:
				fH2 = open ('{}{}/data_clinical_sample.txt'.format(path,Dir),'r')
			lines2 = fH2.readlines()
			tree=0;hL=0
			if D[Dir] == 'Process':
				IDX=3
				for L in range(0,len(lines2),++1):
					tokens = lines2[L].split('\t')
					if tokens[1] == 'SAMPLE_ID' or tokens[1] == 'Sample Identifier' or tokens[1] == '#Sample Identifier': 
						IDX=1
						hL = L+1
						for x in range (0,len(tokens),++1):
							if tokens[x].replace('\n','') == 'ONCOTREE_CODE'or tokens[x].replace('\n','') == 'Oncotree Code':
								tree=x
					elif tokens[0] =='SAMPLE_ID'or tokens[0] == 'Sample Identifier'or tokens[0] == '#Sample Identifier':
						IDX=0
						hL = L+1
						for x in range (0,len(tokens),++1):
							if tokens[x].replace('\n','') == 'ONCOTREE_CODE' or tokens[x].replace('\n','') == 'Oncotree Code':
								tree=x
					elif tokens[2] =='SAMPLE_ID'or tokens[2] == 'Sample Identifier'or tokens[2] == '#Sample Identifier':
						IDX=2
						hL = L+1
						for x in range (0,len(tokens),++1):
							if tokens[x].replace('\n','') == 'ONCOTREE_CODE' or tokens[x].replace('\n','') == 'Oncotree Code':
								tree=x
				for L in range (hL,len(lines2),++1):
					tokens = lines2[L].split('\t')
					if IDX==0:					
						cBio['CSamples'][Dir]={tokens[IDX]:''}
					elif IDX==1:
						cBio['CSamples'][Dir]={tokens[IDX]:''}
					elif IDX==2:
						cBio['CSamples'][Dir]={tokens[IDX]:''}
					else:
						print ('Error:Could not find Sample ID index in {}...'.format(Dir));print (IDX)
						exit(0)
				for L in range (hL,len(lines2),++1):
					tokens = lines2[L].split('\t')
					if IDX==0:					
						cBio['CSamples'][Dir].update({tokens[IDX]:tokens[tree].replace('\n','')})
					elif IDX==1:
						cBio['CSamples'][Dir].update({tokens[IDX]:tokens[tree].replace('\n','')})
					elif IDX==2:
						cBio['CSamples'][Dir].update({tokens[IDX]:tokens[tree].replace('\n','')})
					else:
						print ('Error:Could not find Sample ID index in {}...'.format(Dir));print (IDX)
						exit(0)
				#print ('===============')
			else:#problematic files...
				print ('{}\t{}'.format(Dir,D[Dir]))#;continue 
		except:
			if Dir in ProbFiles:#get filtered files...
				fH2 = open ('{}{}/data_clinical_samples.filtered.txt'.format(path,Dir),'r')#;print ('######{}'.format(Dir))
			else:
				fH2 = open ('{}{}/data_bcr_clinical_data_sample.txt'.format(path,Dir),'r')
			lines2 = fH2.readlines()
			tree=0;hL=0
			if D[Dir] == 'Process':
				IDX=3
				for L in range(0,len(lines2),++1):
					tokens = lines2[L].split('\t')
					if tokens[1] == 'SAMPLE_ID'or tokens[1] == 'Sample Identifier'or tokens[1] == '#Sample Identifier': 
						IDX=1
						hL = L+1
						for x in range (0,len(tokens),++1):
							if tokens[x].replace('\n','') == 'ONCOTREE_CODE' or tokens[x].replace('\n','') == 'Oncotree Code':
								tree=x
					elif tokens[0] =='SAMPLE_ID'or tokens[0] == 'Sample Identifier'or tokens[0] == '#Sample Identifier':
						IDX=0
						hL = L+1
						for x in range (0,len(tokens),++1):
							if tokens[x].replace('\n','') == 'ONCOTREE_CODE'or tokens[x].replace('\n','') == 'Oncotree Code':
								tree=x
					elif tokens[2] =='SAMPLE_ID'or tokens[2] == 'Sample Identifier'or tokens[2] == '#Sample Identifier':
						IDX=2
						hL = L+1
						for x in range (0,len(tokens),++1):
							if tokens[x].replace('\n','') == 'ONCOTREE_CODE'or tokens[x].replace('\n','') == 'Oncotree Code':
								tree=x
				for L in range (hL,len(lines2),++1):
					tokens = lines2[L].split('\t')
					if IDX==0:					
						cBio['CSamples'][Dir]={tokens[IDX]:''}
					elif IDX==1:
						cBio['CSamples'][Dir]={tokens[IDX]:''}
					elif IDX==2:
						cBio['CSamples'][Dir]={tokens[IDX]:''}
					else:
						print ('Error:Could not find Sample ID index in {}...'.format(Dir));print (IDX)
						exit(0)
				for L in range (hL,len(lines2),++1):
					tokens = lines2[L].split('\t')
					if IDX==0:					
						cBio['CSamples'][Dir].update({tokens[IDX]:tokens[tree].replace('\n','')})
					elif IDX==1:
						cBio['CSamples'][Dir].update({tokens[IDX]:tokens[tree].replace('\n','')})
					elif IDX==2:
						cBio['CSamples'][Dir].update({tokens[IDX]:tokens[tree].replace('\n','')})
					else:
						print ('Error:Could not find Sample ID index in {}...'.format(Dir));print (IDX)
						exit(0)
				#print ('===============')
			else:#problematic files...
				print ('{}\t{}'.format(Dir,D[Dir]))#continue
		try:#meta study files...
			fH3 = open ('{}{}/meta_study.txt'.format(path,Dir),'r')
			lines3 = fH3.readlines()
			tempD={}
			for i in info:
				tempD[i]='NA'
			for line in lines3:
				tokens = line.split(':')
				if tokens[0] in info:
					tempD[tokens[0]] = tokens[1].replace('\n','')
				
			
			cBio['MStudy'][Dir]=tempD
				
		except:
			print ('Error:Could not find meta_study.txt in {}'.format(Dir))
			exit(0)
	else:
		print ('Error:Could not find the directory:{}...'.format(Dir))
		exit(0)
#exit(0)
print ('NOTICE:Found {} ready to process files...'.format(len(cBio['CSamples'])))
fOut = open ('cBioPortal.db.tmp.csv','w')
H = '\t'.join(headers)
fOut.write('{}\tgroups\tpmid\tTumorType\n'.format(H))
X=0;OUT=open ('ERROR.txt','w');ED=[]
for f in cBio['MExtended']:
	if f in cBio['CSamples']:
		for sampleID in cBio['MExtended'][f]:
			X += len(cBio['MExtended'][f][sampleID])
			for var in cBio['MExtended'][f][sampleID]:
				tokens = var.split(',')
				line = tokens[4]+'\t'+tokens[0]+'\t'+tokens[6]+'\t'+tokens[1]+'\t'+tokens[8]+'\t'+tokens[9]+'\t'+tokens[5]+'\t'+tokens[7]+'\t'+tokens[2]+'\t'+tokens[3]+'\t'+cBio['MStudy'][f]['groups']+'\t'+cBio['MStudy'][f]['pmid']
				try:
					fOut.write('{}\t{}\n'.format(line,cBio['CSamples'][f][sampleID]))
				except Exception as e:#problematic/low quality samples...
					X-=1
					tmp = str(f)+'\t'+str(e);ED.append(tmp)
	#if f in cBio['MStudy']:
	#	print (f);print(cBio['MStudy'][f]);print('==================')
	else:
		continue
#Reporting problematic/low quality samples...
ED=list(set(ED))
for k in ED:
	OUT.write('{}\n'.format(k))
print ('NOTICE:Done {} variants were written into output file...'.format(X))
#	print (f)
#	print (len( cBio['CSamples'][f]))
#	print (cBio['CSamples'][f])
#	for k in cBio['CSamples'][f]:
#		print (k)
#	exit(0)	





