#/bin/bash

## Software requirments

# python 2.7
# Pysam python module
# Psutil python module
# Biopython
# Perl
# Gmap aligner (>= 2014-12-31)
# R (>= 3.4.0)
# R packages (Sqanti_QC function):
	# ggplot2
	# scales
	# reshape
	# gridExtra
	# grid
# R packages (Sqanti filter function):
	# rpart
	# ROCR
	# caret
	# optparse
	# ggplot2
	# lattice
	# foreach
	# e1071
	# randomForest
	# partykit
	# ipred
	# rpart.plot
	# doMC
	# nnet
	# ROSE
	# pROC
	# MLmetrics


## Setting the working dir

cd example_data_chr16

## Data files

isoformsFasta="all.good.5merge.collapsed.longest_rep16.fasta"
isoExpression="expressionMatrixRSEM16.txt"
fl_pacbio="FL_pacBio"    # path to folder or comma-separeted list of files
sj_covIllumina="SJcoverageIllumina" # path to folder or comma-separeted list of files

# Reference files

refGenome="referenceFiles/mm10_chr16.fa" 
refGTF="referenceFiles/refseq_ensembl16.gtf"
gmapIndex="referenceFiles/mm10_chr16_index/mm10_chr16_index"  # not provided with the example dataset


# Running SQANTI QC 

python ../sqanti_qc.py $isoformsFasta $refGTF $refGenome -fl $fl_pacbio -c $sj_covIllumina -e $isoExpression -n -x $gmapIndex

# Running SQANTI_QC with GTF file (instead of fasta. Not alignment step)

mkdir sqanti_output_fasta
python ../sqanti_qc.py all.good.5merge.collapsed.longest_rep16_corrected.gtf $refGTF $refGenome -fl $fl_pacbio -c $sj_covIllumina -e $isoExpression -n -d sqanti_output_fasta -g


# Running SQANTI filter

SQANTIclass="all.good.5merge.collapsed.longest_rep16_classification.txt"
corrFasta="all.good.5merge.collapsed.longest_rep16_corrected.fasta"

python sqanti_filter.py  $SQANTIclass -i $corrFasta



