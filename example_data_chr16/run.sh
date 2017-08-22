#/bin/bash

## Software requirments

# python 2.7
# Python 2.7
# Pysam python module
# Psutil python module
# Perl
# Gmap aligner (>= 2014-12-31)
# R (>= 3.4.0)

## Setting the working dir

cd exampleData_chr16

## Data files

isoformsFasta="all.good.5merge.collapsed.longest_rep16.fasta"
isoExpression="expressionMatrixRSEM16.txt"
fl_pacbio="FL_pacBio"    # path to folder or comma-separeted list of files
sj_covIllumina="SJcoverageIllumina" # path to folder or comma-separeted list of files

# Reference files

refGenome="referenceFiles/mm10_chr16.fa" 
refGTF="referenceFiles/refseq_ensembl16.gtf"
gmapIndex="referenceFiles/mm10_chr16_index/mm10_chr16_index"  


# Running SQANTI QC 

python ../sqanti_qc.py $isoformsFasta $refGTF $refGenome -fl $fl_pacBio -c $sj_covIllumina -e $isoExpression -n -x $gmapIndex


# Running SQANTI filter

SQANTIclass="all.good.5merge.collapsed.longest_rep16_classification.txt"
corrFasta = "all.good.5merge.collapsed.longest_rep16_corrected.fasta"

python sqanti_filter.py  $SQANTIclass -i $corrFasta



