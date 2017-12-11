# **Welcome to SQANTI: Structural and Quality Annotation of Novel Transcript Isoforms** #


### You can find the preprint of the SQANTI paper at BiorXiv:
* http://biorxiv.org/content/early/2017/03/18/118083 

##**Summary**##

SQANTI is a pipeline for the in-depth characterization of isoforms obtained by full-length transcript sequencing, which are commonly returned in a fasta file format without any extra information about gene/transcript annotation or attribute description. SQANTI provides a wide range of descriptors of transcript quality and generates a graphical report to aid in the interpretation of the sequencing results.

Although SQANTI is oriented to be used for characterization of isoforms generated by PacBio Iso-Seq pipeline, it can be used for any set of transcript isoforms in fasta file format. Besides, it can be applied to any organism.

Moreover, SQANTI adds another functionality, SQANTI filtering function, that allows to filter out artifact transcripts by taking advantage of SQANTI-provided attributes and Machine Learning methods.

SQANTI pipeline steps:

1. First, as long-read sequencing usually has a high rate of errors along sequences, it performs a **reference-based correction of sequences**.

2. Secondly, it **generates genes models** and **classifies transcripts based on splice junctions**.

3. Third, it **predicts ORFs** for each transcript, obtaining information about the coding potential of each sequence. 

4. Finally, it carries out a **deep characterization of isoforms at both transcript and junction level** and **generates a report** with several plots describing in detail the mayor attributes that catalog your set of sequenced isoforms.

5. Together with SQANTI_qc function, the user can use the **SQANTI filter function to remove isoforms potential to be artifacts**. To get this curanted transcriptome, SQANTI filtering uses machine learning methods together with SQANTI_qc attributes to create a classifier of artifacts.


![SQANTI pipeline](sqanti_info/sqanti_pipeline.png?at=master)




## **Required software** ##

* Python 2.7
* Pysam python module
* Psutil python module
* Perl
* Gmap aligner (>= 2014-12-31)
* R (>= 3.4.0)
* R packages (Sqanti_QC function):
        ggplot2,
        scales,
        reshape,
        gridExtra,
        grid.
* R packages (Sqanti_filter function):
        rpart,
        ROCR,
        caret,
        optparse,
        ggplot2,
        lattice,
        foreach,
        e1071,
        randomForest,
        partykit,
        ipred,
        rpart.plot,
        doMC,
        nnet,
        ROSE,
        pROC,
        MLmetrics


## **Running** ##

SQANTI is a program written in python. It is divided into two functions:

* **Sqanti_qc.py**, performing the in-depth characterization of transcripts
* **Sqanti_filter.py**, applying matching learning methods to filter transcripts that are likely to be artifacts.

Below you can see the help page of SQANTI functions where their mandatory and optional arguments are explained.


##**SQANTI_qc.py**##

```
#!bash
usage: sqanti_qc.py [-h] [-g] [-e EXPRESSION] [-x GMAP_INDEX]
                    [-t GMAP_THREADS] [-o OUTPUT] [-d DIR] [-c COVERAGE]
                    [-s SITES] [-w WINDOW] [-n] [-fl FL_COUNT] [-v]
                    isoforms annotation genome

Structural and Quality Annotation of Novel Transcript Isoforms

positional arguments:
  isoforms              Isoforms (Fasta or gtf format; By default "fasta". GTF
                        if specified -g option)
  annotation            Reference annotation file (GTF format)
  genome                Reference genome (Fasta format)

optional arguments:
  -h, --help            show this help message and exit
  -g, --gtf             Use when running Sqanti by using as input a gtf of
                        isoforms
  -e EXPRESSION, --expression EXPRESSION
                        Expression matrix
  -x GMAP_INDEX, --gmap_index GMAP_INDEX
                        Path and prefix of the reference index created by
                        gmap_build. Mandatory unless -g option is specified.
  -t GMAP_THREADS, --gmap_threads GMAP_THREADS
                        Number of threads used during sequence alignment by
                        gmap.
  -o OUTPUT, --output OUTPUT
                        Prefix for output files.
  -d DIR, --dir DIR     Directory for output files. Default: Directory where
                        the script was run.
  -c COVERAGE, --coverage COVERAGE
                        Junction coverage files (comma-separated list of
                        SJ.out.tab files generated by STAR or directory where
                        they are in).
  -s SITES, --sites SITES
                        Set of splice sites to be considered as canonical
                        (comma-separated list of splice sites). Default:
                        GTAG,GCAG,ATAC.
  -w WINDOW, --window WINDOW
                        Size of the window in the genomic DNA screened for
                        Adenine content downstream of TTS
  -n, --name            Use gene_name tag from GTF to define genes. Default:
                        gene_id used to define genes
  -fl FL_COUNT, --fl_count FL_COUNT
                        Full-length PacBio abundance files (comma-separated
                        list of PacBio abundance files generated by PacBio or
                        directory where there are in).
  -v, --version         Display program version number.
```



## Mandatory Input Files: ##

**1. Isoforms:** Isoforms to characterize must be given in a fasta file format. The accepted formats for the fasta headers are:

>\>PB.3.1|chr1:4857776-4897906(+)|i1_c1034/f143p70/2586 

where "PB.3.1" is taken as the isoform ID (typical format returned by PacBio sequencing).

>\>gi|737676268|ref|NM_001545.2| Description 

where NM_001545.2 is taken as the isoform ID (typical format for RefSeq entries).

> \>ENST00000434970.2 cdna chromosome:GRCh38:7:142786213:142786224:1 gene:ENSG00000282431.1 

where ENST00000434970.2 is taken as the isoform ID (typical format for Ensembl entries).

You can also provide Isoforms in gtf file format specifying -g option. "transcript_id" and "gene_id" attributes are mandatory inside GTF file.

# 
**2. Annotation:** To create gene models and therefore annotate each isoform to the corresponding reference gene and reference transcript it's mandatory to provide a reference annotation file in GTF format (http://www.ensembl.org/info/website/upload/gff.html). "transcript_id" and "gene_id" attributes are mandatory for each defined isoform.

# 
**3. Genome:** It must be given in a unique fasta file. "fai" index file must be located in same folder as the fasta file.


#
##Optional Input Files:##

#
**1. Expression Matrix:** Isoform expression. The format must be a tabulated file where the first column corresponds to the isoform identifiers. Next columns must represent the expression levels for each studied sample. The header must have a label for each studied sample.
Expression values can be computed with any software design to calculate transcript expression. We recommend to use short-reads sequencing followed by programs as RSEM, Kallisto or eXpress which allow to calculate accurately isoform expression. 

#
**2. Coverage:** Short-read coverage across junctions represents an informative measure for quality control of sequenced transcripts. STAR aligner output files are the required format files of SQANTI. More specifically, "SJ.out.tab" are the required ones ([STAR manual](http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf)). To get them, short-read data must be align by using STAR aligner. Each sample must be run separately. 

#
**3. FL_count:** As SQANTI was developed for the annotation and quality control of IsoSeq PacBio isoforms, you can provide PacBio abundance files where the number of Full-length reads associated to each PacBio-defined isoform is stored. Specifically, SQANTI uses the "*.abundance.txt" files.

#
**4. Gmap indexes:** Pre-computed genome indexes for GMAP aligner must be provided. 

#

**Note: You can check format of input files in "example_dataset" folder which contains an example.**


##Output Files:##

**1. *name*_corrected.fasta:** Transcript sequences after correction using genome sequence.

**2. *name*_corrected.gff**, **name_corrected.gtf** and **name_corrected.sam**: Alignment of corrected sequences 

**3. *name*.faa:** File with predicted ORF sequences for corrected isoforms. GeneMarkS-T (http://exon.gatech.edu/GeneMark/) was used to predict coding sequences in eukaryotic transcripts. Also a GMST folder is created with extra information returned by GMST.

**4. RTS folder:** Folder which stores a file with information about junctions that are potential to have suffered RT-switching during retrotranscription. RT switching algorithm locates 8 nt direct repeats characteristic of RT switching between the upstream mRNA boundary of the non-canonical intron and the intron region adjacent to the downstream exon boundary with a slide of 1 nucleotide upstream/downstream. FSM transcripts with the highest mean expression in each gene are assumed to serve as templates for RT switching and are excluded from the analysis.

**5. *name*_classification.txt:** File with attribute information at isoform level (table explaining feature meaning inside output_info).

**6. *name*_junctions.txt:** File with attribute information at splice junction level (table explaining feature meaning inside output_info).

**7. *name*_report.pdf:** PDF file showing different quality control and descriptive plots. An example can be found [here:](https://bitbucket.org/ConesaLab/sqanti/src/a0f3d2a16452f304d91c3fb0020d979284cb82b3/example_data_chr16/all.good.5merge.collapsed.longest_rep16_report.pdf?at=master)



##**SQANTI_filter.py**##

```
#!bash
usage: sqanti_filter.py [-h] [-t T] [-a A] [-i FASTA] [-d DIR]
                        [-p TP_ISOFORMS] [-n TN_ISOFORMS] [-v]
                        sqanti_class

Filtering of Isoforms based on SQANTI attributes

positional arguments:
  sqanti_class          SQANTI classification at isoform level

optional arguments:
  -h, --help            show this help message and exit
  -t T, -ml_threshold T
                        Machine learning probability threshold to classify
                        posiive isoforms
  -a A, -intrapriming A
                        Adenine percentage at genomic 3' end to flag an
                        isoform as intra-priming
  -i FASTA, --isoforms FASTA
                        Fasta isoform file to be filtered by SQANTI
  -d DIR, --dir DIR     Output directory name. Default: "Filter_out" directory
                        where the script was run
  -p TP_ISOFORMS, -TP_isoforms TP_ISOFORMS
                        File with true positive set of isoforms
  -n TN_ISOFORMS, -TN_isoforms TN_ISOFORMS
                        File with true negative set of isoforms
  -v, --version         Display program version number.

```
## Mandatory Input Files: ##

**1. SQANTI classification file:** SQANTI output file with attribute information at isoform level.

##Optional Input Files:##

**1. Corrected isoforms file:** Isoforms to be filtered can be provided in a fasta file format to be automatically filtered by sqanti filtering.


##Output Files:##

**1. Confusion_matrix_trainingset.txt:** confusion matrix obtained on the cross-validation (10 fold, repeated 10 times) preformed on the training set. The rows correspond to the predictions, the columns corresponds to the reference (TP_isoforms, TN_isoforms).

**2. Variable_Importance_trainingset.txt:** Importance of each variable on the final classifier. Values obtained with the function varImp of the caret library.

**3. ROC_curves_RF.pdf:** ROC curves obtained on the cross-validation. One doted ROC curve by repetition. The red curve corresponds to the curve computed from the 10 repetitions. The area under ROC is computed and written at the bottom right of the plot.

**4. feature_importance_training.pdf:** boxplot on the training set of the values of the 5 variables used by the classifier in function of the class (TP/TN). Log scale. In order to avoid infinite values in log scale we added plus one to each value for the variables “bite”, “minimum coverage”, “minimum sample coverage”, “isoform expression”.

**5. feature_importance_ total.pdf:** boxplot of the values of the 5 variables used by the classifier on the predictions performed on the multi-exonic transcripts. Log scale. In order to avoid infinite value in log scale, we added plus one to each value for the variables “bite”, “minimum coverage”, “minimum sample coverage”, “isoform expression”.

**6. curated_transcriptome.txt:**  One columns file providing the set of isoforms don't predicted as artifacs.


##**Version**##

* SQANTI 1.2


##**Contact Information**##

* Lorena de la Fuente.
lfuente@cipf.es

* Manuel Tardaguila.
manueltar@ufl.edu




