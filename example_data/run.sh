#/bin/bash

##### To run SQANTI you need to download the mm10.fa genome and generate the mm10_index for gmap alignment. Not incorporated into example folder because of size.

# Setting the working dir

cd SQANTI_exampleData

# Data files

isoformsFasta="all.good.5merge.collapsed.longest_rep.fasta"
isoExpression="expressionMatrixRSEM.txt"
fl_pacbio="FL_pacbio"    # path to folder or comma-separeted list of files
sj_covIllumina="SJcoverageIllumina" # path to folder or comma-separeted list of files
tp_iso="TP_TN_iso/TP.txt"
tn_iso="TP_TN_iso/TN.txt"


# Reference files

refGenome="referenceFiles/mm10.fa" # Not in downloadable example dataset 
refGTF="referenceFiles/refseq_ensembl_smallRNAsfiltered.gtf"
gmapIndex="referenceFiles/mm10_index/mm10_index"  # Not in downloadable example dataset 


# Running SQANTI with correction step

python sqanti.py  $isoformsFasta $refGTF $isoExpression $refGenome -fl $fl_pacbio -c $sj_covIllumina -n -x $gmapIndex

# If SQANTI was previously run, we can use the *corrected* files generated during the previous run to get attributes (option -m). Important if we just want to change expression/coverage data but using the same transcriptome. Not need of specified the gmap index.

python sqanti.py  $isoformsFasta $refGTF $isoExpression $refGenome -fl $fl_pacbio -c $sj_covIllumina -n -m

# Running the SQANTI filter function to create the classifier of isoforms

SQANTIclass_output="all.good.5merge.collapsed.longest_rep_classification.txt"
filter_outputDir="myclassifier_test"

python sqanti_filtering.py $tp_iso $tn_iso $SQANTIclass_output -i $isoformsFasta -d $filter_outputDir



