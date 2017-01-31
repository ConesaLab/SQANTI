# **Welcome to SQANTI: Structural and Quality Annotation of Novel Transcript Isoforms** #

 
##**Summary**##


SQANTI is a pipeline for the in-depth characterization of isoforms obtained by full-length transcript sequencing, which are commonly returned in a fasta file format without any extra information about gene/transcript annotation or attribute description. SQANTI provides a wide range of descriptors of transcript quality and generates a graphical report to aid in the interpretation of the sequencing results.

Although SQANTI is oriented to be used for characterization of isoforms generated by PacBio Iso-Seq pipeline, it can be used for any set of transcript isoforms in fasta file format. Besides it can be applied to any organism.

SQANTI pipeline steps:

1. First, as long-read sequencing usually has a high rate of errors along sequences, it performs a **reference-based correction of sequences**.

2. Secondly, **generates genes models** and **classifies transcripts based on splice junctions**.
      ![Transcript Classification](https://bitbucket.org/repo/kpnA5g/images/803890880-Diapositiva1.png)

3. Third, **computes a ORF prediction** for each isoform, obtaining information about the coding potential of each sequenced isoform. 

4. Finally it carries out a **deep characterization of isoforms at both transcript and junction level** and **generates a report** with several plots describing in detail the mayor attributes that catalog your set of sequenced-isoforms.

![SQANTI WORKFLOW](https://bitbucket.org/repo/kpnA5g/images/3187661279-Diapositiva8.png)

## **Required software** ##

* Python
* Pysam python module
* Perl
* Gmap aligner (version....)
* R environment


## **Running** ##

SQANTI is a program written in python. To run it use *sqanti.py*. 

Below you can see the help page of SQANTI where its mandatory and optional arguments are explained:

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


##**Input files**##


### Mandatory Input Files ###

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
**4. Expression Matrix:** SQANTI needs information about the pre-computed expression of each given isoform to calculate some quality control attributes. The format must be tabulated file where the first column corresponds to the isoform identifiers found in the fasta file. Next columns must represent the expression levels for each studied sample. Expression values can be computed with any software design to calculate transcript expression. We recommend to use short-reads sequencing followed by programs as RSEM, Kallisto or eXpress which allow to calculate accurately isoform expression. 

#
### Optional Input Files ###

**1. Coverage:** Short-read coverage in junctions represents an informative measure for quality control of sequenced transcripts. STAR aligner output files are the required format files of SQANTI. More specifically, "SJ.out.tab" are the required ones ([STAR manual](http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf)). To get them, short-read data must be align by using STAR aligner. Each sample must be run separately. 

#
**2. FL_count:** As SQANTI was developed for the annotation and quality control of IsoSeq PacBio isoforms, you can provide PacBio abundance files where the number of Full-length reads associated to each PacBio-defined isoform is stored. Specifically, SQANTI uses the "*.abundance.txt" files.

#
**3. Gmap indexes: --------
#
*Note: You can see file format examples in "example" folder.*


##**Output files**##

##**Version**##

* SQANTI 0.1


##**Contact Information**##

Lorena de la Fuente 
lfuente@cipf.es

Genomics of Gene Expression Lab
Centro de Investigación Príncipe Felipe
Eduardo Primo Yúfera, 3 46012 Valencia (Spain)