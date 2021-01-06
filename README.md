# CACSV
Clinically Actionable Cancer Somatic Variants (CACSV) is genetics dataset that is built of publicly available cancer mutations simulated from different tumor sites and catalogued based on the AMP-ASCO-CAP 2017 recommendations. 
## PREREQUISITE
--cancer databases--
```
1. You need to install Python >= 2.7
2. You need to download the latest version of the following databases:
(a) COSMIC
CosmicMutantExportCensus.tsv  (https://cancer.sanger.ac.uk/cosmic/download )
cancer_gene_census.csv  (https://cancer.sanger.ac.uk/cosmic/download )
(b) cBioPortal
data_mutations_extended.txt files per study  (https://www.cbioportal.org/datasets)
(c) OncoKB
allActionableVariants.txt  (https://www.oncokb.org/apiAccess )
(d) CCGD
CCGD_export.csv (http://ccgd-starrlab.oit.umn.edu/search.html )
```

## USAGE
--prepare the files--
```
python PrepCCGD.py CCGD_export.csv
python CreatMasterFile.py cancer_gene_census.csv,CCGD_export.db.csv TypesMatching.csv,TypesMatching.CCGD.csv
python Collection.2.py TypesMatching_cBioPortal.csv cBioPortal_all/datahub/public/[full path]
python PrePcBio.2.py cBioPortal.db.1.csv,Tumor_Type.3.csv
python SV-Collection.COSMIC.V4.2.py CosmicMutantExportCensus.mmddyyy.tsv Tumor_Type.3.csv
python mAke.oncoKB_based_db.1.py allActionableVariants.txt
--move the generated files to final_db directory--

```
--run the classifier--
```
python Classifier.3.5.py somatic_variants_text_file 
 
```
## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
