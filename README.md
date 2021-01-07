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
python3 CreatMasterFile.py cancer_gene_census.csv,CCGD_export.db.csv TypesMatching.csv,TypesMatching.CCGD.csv
python3 Collection.2.py TypesMatching_cBioPortal.csv cBioPortal_all/datahub/public/[full path]
python3 PrePcBio.2.py cBioPortal.db.1.csv,Tumor_Type.3.csv
python3 SV-Collection.COSMIC.V4.2.py CosmicMutantExportCensus.mmddyyy.tsv Tumor_Type.3.csv
python3 mAke.oncoKB_based_db.1.py allActionableVariants.txt
--move the generated files to final_db directory--

```
--run the classifier--
```
python3 Classifier.3.5.py somatic_variants_text_file 
 
```
## LICENSE FOR DATASET

CACSV dataset is free for non-commercial use without warranty. Please contact the authors for commercial use.

## LICENSE ON CODE

This license only concerns the code fully written by us.

The MIT License (MIT)

Copyright (c) 2021, Turki M. Sobahy (KFSHRC-J)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
