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
      ![Transcript Classification](https://bitbucket.org/repo/kpnA5g/images/803890880-Diapositiva1.png)

3. Third, it **predicts ORFs** for each transcript, obtaining information about the coding potential of each sequence. 

4. Finally, it carries out a **deep characterization of isoforms at both transcript and junction level** and **generates a report** with several plots describing in detail the mayor attributes that catalog your set of sequenced isoforms.

5. Together with SQANTI function, the user can use the **SQANTI filtering function to remove isoforms potential to be artifacts**. To get this curanted transcriptome, SQANTI filtering uses machine learning methods together with SQANTI-defined attributes to create a classifier of artifacts.


![SQANTI WORKFLOW](https://bitbucket.org/repo/kpnA5g/images/2483959904-Diapositiva8.png)

## **Required software** ##

* Python 2.7
* Pysam python module
* Psutil python module
* Perl
* Gmap aligner (>= 2014-12-31)
* R (>= 3.3.2)


## **Running** ##

SQANTI is a program written in python. It is divided into two functions:

* **Sqanti.py**, performing the in-depth characterization of transcripts
* **Sqanti_filtering.py**, applying matching learning methods to filter transcripts that are likely to be artifacts.

Below you can see the help page of SQANTI functions where their mandatory and optional arguments are explained.


##**SQANTI.py**##

```
#!bash
usage: sqanti.py [-h] [-x GMAP_INDEX] [-o OUTPUT] [-d DIR] [-c COVERAGE]
                 [-s SITES] [-n] [-fl FL_COUNT] [-m] [-v]
                 isoforms annotation expression genome

Structural and Quality Annotation of Novel Transcript Isoforms

positional arguments:
  isoforms              Isoforms (Fasta format)
  annotation            Reference annotation file (GTF format)
  expression            Expression matrix
  genome                Reference genome (Fasta format)

optional arguments:
  -h, --help            show this help message and exit
  -x GMAP_INDEX, --gmap_index GMAP_INDEX
                        Path and prefix of the reference index created by
                        gmap_build. Mandatory unless -m option is specified.
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
  -n, --name            Use gene_name tag from GTF to define genes. Default:
                        gene_id used to define genes
  -fl FL_COUNT, --fl_count FL_COUNT
                        Full-length PacBio abundance files (comma-separated
                        list of PacBio abundance files generated by PacBio or
                        directory where there are in).
  -m, --mode            Use to run Sqanti when gtf and faa files have been
                        already created in previous runs.
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

# 
**2. Annotation:** To create gene models and therefore annotate each isoform to the corresponding reference gene and reference transcript it's mandatory to provide a reference annotation file in GTF format (http://www.ensembl.org/info/website/upload/gff.html).

*Note: It's recommended to discard smallRNAs (miRNAs, snoRNA, etc) to avoid a misannotation of novel genes.*
# 
**3. Genome:** Mandatory for the error-correction of sequences during genome mapping. It must be given in a unique fasta file.

#
**4. Expression Matrix:** SQANTI needs information about the pre-computed expression of each given isoform to calculate some quality control attributes. The format must be tabulated file where the first column corresponds to the isoform identifiers. Next columns must represent the expression levels for each studied sample. The header must have a label for each studied sample.
Expression values can be computed with any software design to calculate transcript expression. We recommend to use short-reads sequencing followed by programs as RSEM, Kallisto or eXpress which allow to calculate accurately isoform expression. 

#
##Optional Input Files:##

**1. Coverage:** Short-read coverage across junctions represents an informative measure for quality control of sequenced transcripts. STAR aligner output files are the required format files of SQANTI. More specifically, "SJ.out.tab" are the required ones ([STAR manual](http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf)). To get them, short-read data must be align by using STAR aligner. Each sample must be run separately. 

#
**2. FL_count:** As SQANTI was developed for the annotation and quality control of IsoSeq PacBio isoforms, you can provide PacBio abundance files where the number of Full-length reads associated to each PacBio-defined isoform is stored. Specifically, SQANTI uses the "*.abundance.txt" files.

#
**3. Gmap indexes:** Pre-computed genome indexes for GMAP aligner must be provided. 

#

**Note: You can check format of input files in "example_dataset" folder which contains a real example.**


##Output Files:##

**1. *name*_corrected.fasta:** Transcript sequences after correction using genome sequence.

**2. *name*_corrected.gff**, **name_corrected.gtf** and **name_corrected.sam**: Alignment of corrected sequences 

**3. *name*.faa:** File with predicted ORF sequences for corrected isoforms. GeneMarkS-T (http://exon.gatech.edu/GeneMark/) was used to predict coding sequences in eukaryotic transcripts. Also a GMST folder is created with extra information returned by GMST.

**4. RTS folder:** Folder which stores a file with information about junctions that are potential to have suffered RT-switching during retrotranscription. RT switching algorithm locates 8 nt direct repeats characteristic of RT switching between the upstream mRNA boundary of the non-canonical intron and the intron region adjacent to the downstream exon boundary with a slide of 1 nucleotide upstream/downstream. FSM transcripts with the highest mean expression in each gene are assumed to serve as templates for RT switching and are excluded from the analysis.

**5. *name*_classification.txt:** File with attribute information at isoform level (table explaining features at the bottom of the file)

**6. *name*_junctions.txt:** File with attribute information at splice junction level.

**7. *name*_Report.pdf:** PDF file showing different quality control and descriptive plots. 



##**SQANTI_filtering.py**##

```
#!bash
usage: sqanti_filtering.py [-h] [-i FASTA] [-d DIR] [-v]
                           TP_isoforms TN_isoforms sqanti_class

Filtering of Isoforms based on SQANTI attributes

positional arguments:
  TP_isoforms           File with true positive set of isoforms
  TN_isoforms           File with true negative set of isoforms
  sqanti_class          SQANTI classification at isoform level

optional arguments:
  -h, --help            show this help message and exit
  -i FASTA, --isoforms FASTA
                        Fasta isoform file to be filtered by SQANTI filtering
  -d DIR, --dir DIR     Output directory name. Default: "Classifier_out"
                        directory where the script was run.
  -v, --version         Display program version number.

```
## Mandatory Input Files: ##

**1. TP_isoforms:** One column file with true positive set of isoforms.

**2. TN_isoforms:** One column file with true negative set of isoforms.

**3. SQANTI classification file:** SQANTI output file with attribute information at isoform level.

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

* SQANTI 0.1


##**Contact Information**##

* Lorena de la Fuente.
lfuente@cipf.es

* Manuel Tardaguila.
manueltar@ufl.edu


##**Appendix**##

* Attributes at isoform level defined by SQANTI in **name_classification.txt** file

![table_supp1_HOJA1-1.png](https://bitbucket.org/repo/kpnA5g/images/3507259918-table_supp1_HOJA1-1.png)


* Attributes at splice junction level defined by SQANTI in **name_junctions.txt** file

![table_supp1_HOJA2-1 2.png](https://bitbucket.org/repo/kpnA5g/images/2891116287-table_supp1_HOJA2-1%202.png)