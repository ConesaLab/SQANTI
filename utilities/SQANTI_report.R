#####################################
##### SQANTI report generation ######
#####################################

### Author: Lorena de la Fuente
### Date: 14/10/2016

#********************** Taking arguments from python script

args <- commandArgs(trailingOnly = TRUE)
class.file = args[1]
junc.file = args[2]
report.file = paste(strsplit(class.file, "_classification.txt")[[1]][1], "report.pdf", sep="_")


#********************** Packages (install if not found)

list.of.packages <- c("ggplot2", "scales", "reshape", "gridExtra", "grid")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


library(ggplot2)
library(scales)
library(reshape)
library(gridExtra)
library(grid)


#********************** Reading Information

########## Classification information

sqantiData = read.table(file=class.file, header=T, as.is=T, sep="\t")
rownames(sqantiData) = sqantiData$isoform

if (!all(is.na(sqantiData$iso_exp))){
  sorted <- sqantiData[order(sqantiData$iso_exp, decreasing = T),]
  FSMhighestExpIsoPerGene <- sorted[(!duplicated(sorted$associated_gene) & sorted$structural_category=="full-splice_match"),"isoform"]
  sqantiData[which(sqantiData$isoform%in%FSMhighestExpIsoPerGene),"RTS_stage"] <- FALSE
  write.table(sqantiData, file=class.file, row.names=FALSE, quote=F, sep="\t")
}

xaxislabelsF1 = c("FSM", "ISM", "NIC", "NNC", "Genic\nGenomic",  "Antisense", "Fusion","Intergenic", "Genic\nIntron")
sqantiData$structural_category = factor(sqantiData$structural_category, 
                                       labels = xaxislabelsF1, 
                                       levels = c("full-splice_match","incomplete-splice_match","novel_in_catalog","novel_not_in_catalog", "genic","antisense","fusion","intergenic","genic_intron"), 
                                       ordered=TRUE)


########### Junction information

junctionsData = read.table(file=junc.file, header=T, as.is=T, sep="\t")

if (!all(is.na(sqantiData$iso_exp))){
  junctionsData[which(junctionsData$isoform%in%FSMhighestExpIsoPerGene),"RTS_junction"] <- FALSE
  write.table(junctionsData, file=junc.file, row.names=FALSE, quote=F, sep="\t")
}


########## Generating plots


#*** Global plot parameters

myPalette = c("#6BAED6","#FC8D59","#78C679","coral2","#969696","#66C2A4", "goldenrod1", "darksalmon", "#41B6C4","tomato3", "#FE9929")

mytheme <- theme_classic(base_family = "Helvetica") +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4)) +
  theme(axis.title.x = element_text(size=14, margin=margin(5,0,0,0)),
        axis.text.x  = element_text(margin=margin(7,0,0,0), size=13),
        axis.title.y = element_text(size=14,  margin=margin(0,15,0,0)),
        axis.text.y  = element_text(vjust=0.5, size=13) ) +
  theme(legend.text = element_text(size = 10), legend.title = element_text(size=11), legend.key.size = unit(0.5, "cm")) +
  theme(plot.title = element_text(lineheight=.4, size=13, margin=margin(-20,0,0,0))) +
  #theme(legend.justification=c(1,1), legend.position=c(1,1)) +
  #theme(plot.margin = unit(c(1.5,0.5,1,1), "cm")) +
  theme(plot.margin = unit(c(2.5,1,1,1), "cm")) 


mytheme_bw <- theme_bw(base_family = "Helvetica") +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4)) +
  theme(axis.title.x = element_text(size=14, margin=margin(5,0,0,0)),
        axis.text.x  = element_text(margin=margin(7,0,0,0), size=13),
        axis.title.y = element_text(size=14,  margin=margin(0,15,0,0)),
        axis.text.y  = element_text(vjust=0.5, size=13) ) +
  theme(legend.text = element_text(size = 10), legend.title = element_text(size=11), legend.key.size = unit(0.5, "cm")) +
  theme(plot.title = element_text(lineheight=.4, size=13, margin=margin(-20,0,0,0))) +
  #theme(legend.justification=c(1,1), legend.position=c(1,1)) +
  theme(plot.margin = unit(c(2.5,1,1,1), "cm")) 

roundUp = function(x) 10^ceiling(log10(x))
RoundUp = function(from,to) ceiling(from/to)*to

#*** PLOT 0: Number Isoforms per Gene

sqantiData[grep("novelGene",  sqantiData$associated_gene), "novelGene"] <- "Novel Genes"
sqantiData[-grep("novelGene",  sqantiData$associated_gene), "novelGene"] <- "Annotated Genes"
sqantiData$novelGene = factor(sqantiData$novelGene, 
                                          levels = c("Novel Genes","Annotated Genes"), 
                                          ordered=TRUE)
sqantiData[which(sqantiData$exons>1), "exonCat"] <- "MultiExon"
sqantiData[which(sqantiData$exons==1), "exonCat"] <- "MonoExon"
sqantiData$exonCat = factor(sqantiData$exonCat, 
                                    levels = c("MultiExon","MonoExon"), 
                                    ordered=TRUE)


sqantiData$gene_exp = as.character(sqantiData$gene_exp)


if (!all(is.na(sqantiData$gene_exp))){
  isoPerGene = aggregate(sqantiData$isoform,
                         by = list("associatedGene" = sqantiData$associated_gene,
                                   "novelGene" = sqantiData$novelGene,
                                   "FSM_class" = sqantiData$FSM_class,
                                   "geneExp"=sqantiData$gene_exp),
                         length)
}else{isoPerGene = aggregate(sqantiData$isoform, 
                             by = list("associatedGene" = sqantiData$associated_gene, 
                                       "novelGene" = sqantiData$novelGene, 
                                       "FSM_class" = sqantiData$FSM_class), 
                             length)
}

colnames(isoPerGene)[ncol(isoPerGene)] <- "nIso"


isoPerGene$FSM_class2 = factor(isoPerGene$FSM_class, 
                               levels = c("A", "B", "C"), 
                               labels = c("MonoIsoform Gene", "MultiIsoform Genes\nwithout expression\nof a FSM", "MultiIsoform Genes\nexpressing at least\none FSM"), 
                               ordered=TRUE)

isoPerGene$novelGene = factor(isoPerGene$novelGene, 
                              levels = c("Annotated Genes", "Novel Genes"), 
                              ordered=TRUE)

maxiso = max(isoPerGene$nIso)+1

isoPerGene$range =cut(isoPerGene$nIso, breaks = c(0,1,3,5,maxiso), labels = c("1", "2-3", "4-5", ">=6"))

legendLabelF1 = tools::toTitleCase(gsub("_", " ", levels(as.factor(sqantiData$coding))))

p0 <- ggplot(isoPerGene, aes(x=range, fill=range)) +
  geom_bar(aes(stat = "count", y= (..count..)/sum(..count..)), color="black", size=0.3, width=0.7) +
  guides(fill=FALSE) +
  scale_y_continuous(labels = percent, expand = c(0,0)) +
  scale_fill_manual(values = myPalette[c(2:5)]) +
  labs(x ="# Isoforms per Gene", title="Distribution of isoforms per gene\n\n\n", y = "% Genes") +
  mytheme



#*** PLOT 1: Structural Classification

p1 <- ggplot(data=sqantiData, aes(x=structural_category)) +
  geom_bar(aes(y = (..count..)/sum(..count..), alpha=coding, fill=structural_category), color="black", size=0.3, width=0.7) +
  scale_y_continuous(labels = percent, expand = c(0,0), limits = c(0,1)) +
  geom_text(aes(y = ((..count..)/sum(..count..)), label = scales::percent((..count..)/sum(..count..))), stat = "count", vjust = -0.25)  +
  scale_x_discrete(drop=FALSE) +
  scale_alpha_manual(values=c(1,0.3), 
                     name = "Coding prediction", 
                     labels = legendLabelF1)+
  scale_fill_manual(values = myPalette, guide='none') + 
  xlab("") + 
  ylab("% Transcripts") +
  mytheme + 
  geom_blank(aes(y=((..count..)/sum(..count..))), stat = "count") +
  theme(axis.text.x = element_text(angle = 45)) +
  ggtitle("Isoform distribution across structural categories\n\n" ) +
  theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12)) +
  theme(legend.justification=c(1,1), legend.position=c(1,1)) 


#*** PLOTS 2-3: refLength and refExons for ISM and FSM transcripts. Plot if any ISM or FSM transcript


if (nrow(sqantiData[sqantiData$structural_category%in%c("FSM","ISM"),])!=0){
  
  sqantiFSMISM = sqantiData[sqantiData$structural_category%in%c("FSM","ISM"),]
  sqantiFSMISM$structural_category = factor(sqantiFSMISM$structural_category, 
                                          levels = c("FSM","ISM"), 
                                          ordered=TRUE)

  p2 <- ggplot(data=sqantiFSMISM, aes(x=structural_category, y=ref_length, fill=structural_category)) +
    geom_boxplot(color="black", size=0.3, outlier.size = 0.2) + mytheme +
    scale_fill_manual(values = myPalette) +
    scale_x_discrete(drop=FALSE) +
    guides(fill=FALSE) +
    xlab("") +  
    ylab("Matched Reference Length (bp)") +
    labs(title="Length distribution of matched reference transcripts\n\n\n",  
         subtitle="Just applicable to FSM and ISM categories\n\n") 


  p3 <- ggplot(data=sqantiFSMISM, aes(x=structural_category, y=ref_exons, fill=structural_category)) +
    geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +   
    scale_x_discrete(drop=FALSE) +
    xlab("") +  
    ylab("Matched Reference exon number") +
    scale_fill_manual(values = myPalette) +
    guides(fill=FALSE) +
    mytheme +
    labs(title="Exon number distribution of matched reference transcripts\n\n\n",  
         subtitle="Just applicable to FSM and ISM categories\n\n") 

}


# PLOT 4: Transcript length by category

p4 <- ggplot(data=sqantiData, aes(x=structural_category, y=length, fill=structural_category)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  scale_x_discrete(drop=FALSE) +
  ylab("Transcript Length (bp)") +
  scale_fill_manual(values = myPalette) +
  guides(fill=FALSE) +
  mytheme  + theme(axis.text.x = element_text(angle = 45)) + 
  theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12))+
  ggtitle("Transcript length distribution by structural classification\n\n" ) +  
  theme(axis.title.x=element_blank()) 


# PLOT 5: Exon number by category 

p5 <- ggplot(data=sqantiData, aes(x=structural_category, y=exons, fill=structural_category)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  ylab("Number of exons") +
  scale_x_discrete(drop=FALSE) +
  scale_fill_manual(values = myPalette) +
  guides(fill=FALSE) +
  mytheme  + theme(axis.text.x = element_text(angle = 45)) + 
  theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12))+
  ggtitle("Exon number distribution by structural classification\n\n" ) +
  theme(axis.title.x=element_blank())


# PLOT 6: Exon number vs gene type

p6 <- ggplot(data=sqantiData, aes(x=novelGene)) +
  geom_bar(position="fill",aes(y = (..count..)/sum(..count..), fill=exonCat), color="black", size=0.3, width=0.5) +
  #geom_text(aes(y = (..count..), label=(..count..)), stat = "count", vjust = -0.25)  +
  scale_x_discrete(drop=FALSE) +
  scale_y_continuous(labels = percent, expand = c(0,0)) +
  scale_fill_manual(name = "Transcript type", 
                    values = myPalette[c(2:5)]) +
  ylab("% Transcripts ") +  
  mytheme +
  theme(axis.title.x=element_blank()) + 
  theme(legend.position="bottom") +   
  ggtitle("Distribution of mono/multi exon transcripts\n\n" ) 



# PLOT 7: number of Iso per gene vs gene type

p7 <- ggplot(data=isoPerGene, aes(x=novelGene)) +
  geom_bar(position="fill", aes(y = (..count..)/sum(..count..), fill=range), color="black", size=0.3, width=0.5) +
  scale_y_continuous(labels = percent, expand = c(0,0)) +
  scale_x_discrete(drop=FALSE) +
  scale_fill_manual(name = "# Isoforms Per Gene",
                    values = myPalette[c(2:5)]) +
  ylab("% Genes ") +  
  xlab("Gene Type") +
  mytheme +
  theme(axis.title.x=element_blank()) + 
  theme(legend.position="bottom") +   
  guides(fill = guide_legend(keywidth = 0.9, keyheight = 0.9)) +
  ggtitle("Distribution of number of isoforms\n\n\n\n" ) 


# PLOT 8: Expression, if isoform expression provided

if (!all(is.na(sqantiData$iso_exp))){
  
  p8 <- ggplot(data=sqantiData, aes(x=structural_category, y=log2(iso_exp+1), fill=structural_category)) +
    geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
    scale_x_discrete(drop=FALSE) +
    ylab("log2(transcript expression +1)") +
    scale_fill_manual(values = myPalette) +
    guides(fill=FALSE) +
    mytheme  + theme(axis.text.x = element_text(angle = 45)) + 
    theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12))+
    theme(axis.title.x=element_blank()) + 
    #theme(plot.margin = unit(c(1.5,1,0,1), "cm")) +
    ggtitle("Transcript expression by structural category\n\n" )

}



# PLOT 9: FL number, if FL count provided

if (!all(is.na(sqantiData$FL))){
  
  p9 <- ggplot(data=sqantiData, aes(x=structural_category, y=log2(FL+1), fill=structural_category)) +
    geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
    ylab("log2( # FL reads + 1)") +
    scale_x_discrete(drop=FALSE) +
    scale_fill_manual(values = myPalette) +
    guides(fill=FALSE) +
    mytheme +
    theme(axis.text.x = element_text(angle = 45)) + 
    theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12))+
    theme(axis.title.x=element_blank()) + 
    #theme(plot.margin = unit(c(1.5,1,0,1), "cm")) +
    ggtitle("Number of FL reads per transcript by structural category\n\n" )
  
}


# PLOT 10: Gene Expression, if expresion provided

if (!all(is.na(sqantiData$iso_exp))){

  p10 <- ggplot(data=isoPerGene, aes(x=novelGene, y=log2(geneExp+1), fill=novelGene)) +
    geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
    scale_x_discrete(drop=FALSE) +
    xlab("Structural Classification") +  
    ylab("log2( # Short reads + 1)") +
    scale_fill_manual(values = myPalette[c(3:4)]) +
    guides(fill=FALSE) +
    mytheme +
    theme(axis.title.x=element_blank()) + 
    #theme(plot.margin = unit(c(1.5,1,0.5,1), "cm")) +
    ggtitle("Gene expression by type of gene annotation\n\n" )
}



# PLOT 11: Gene FL number, if FL count provided

if (!all(is.na(sqantiData$FL))){
  
  FL_gene = aggregate(as.integer(sqantiData$FL), by = list("associatedGene" = sqantiData$associated_gene), sum)
  colnames(FL_gene)[ncol(FL_gene)] <- "FL_gene"
  isoPerGene = merge(isoPerGene, FL_gene, by="associatedGene")
  
  p11 <- ggplot(data=isoPerGene, aes(x=novelGene, y=log2(FL_gene+1), fill=novelGene)) +
    geom_boxplot(color="black", size=0.3,outlier.size = 0.2) +
    scale_x_discrete(drop=FALSE) +
    ylab("log2( # FL reads + 1)") +
    scale_fill_manual(values = myPalette[c(3:4)]) +
    guides(fill=FALSE) +
    mytheme +
    theme(axis.title.x=element_blank()) + 
    #theme(plot.margin = unit(c(1.5,1,0.5,1), "cm")) +
    ggtitle("Number of FL reads per Gene by type of gene annotation\n\n" )
  
}



# PLOT 12: NNC expression genes vs not NNC expression genes

# NNC expression genes vs not NNC expression genes


if (!all(is.na(sqantiData$gene_exp))){
  
  if (nrow(sqantiData[sqantiData$structural_category=="NNC",])!=0){
    
    NNC_genes = unique(sqantiData[sqantiData$structural_category=="NNC","associated_gene"]) 
    notNNC_genes = unique(sqantiData[!sqantiData$associated_gene%in%NNC_genes,"associated_gene"]) 
    isoPerGene[isoPerGene$associatedGene %in% notNNC_genes, "NNC_class"] <- "Genes not expressing\n NNC isoforms"
    isoPerGene[isoPerGene$associatedGene %in% NNC_genes, "NNC_class"] <- "Genes expressing\n NNC isoforms"
    
    isoPerGene$NNC_class = factor(isoPerGene$NNC_class, levels=c("Genes expressing\n NNC isoforms","Genes not expressing\n NNC isoforms"),
                            labels=c("Genes expressing\n NNC isoforms","Genes not expressing\n NNC isoforms"), order=T)
    
    p12 <- ggplot(data=isoPerGene[!is.na(isoPerGene$NNC_class),], aes(x=NNC_class, y=log2(geneExp+1), fill=NNC_class)) +
      geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
      xlab("") +  
      ylab("log2(# Short reads per gene + 1)") +
      scale_x_discrete(drop=FALSE) +
      scale_fill_manual(values = c(myPalette[4],"grey38")) +
      guides(fill=FALSE) +
      mytheme +
      theme(axis.title.x=element_blank()) + 
      ggtitle("Gene expression levels between NNC and not NNC containing genes\n\n" ) 
    }
}


# PLOT 13: Genes expression to only FSM Genes, only NNC Genes and both containing genes

if (!all(is.na(sqantiData$gene_exp))){
  
  
  if (nrow(sqantiData[sqantiData$structural_category=="NNC",])!=0 & nrow(sqantiData[sqantiData$structural_category=="FSM",])!=0 ){
    
    FSM_just_genes = unique(sqantiData[sqantiData$FSM_class=="A" & sqantiData$structural_category=="FSM","associated_gene"])
    NNC_just_genes = unique(sqantiData[sqantiData$FSM_class=="A" & sqantiData$structural_category=="NNC","associated_gene"])
    FSMandNNCgenes = unique(sqantiData[sqantiData$FSM_class=="C" & sqantiData$structural_category=="NNC","associated_gene"])
    isoPerGene[isoPerGene$associatedGene %in% FSMandNNCgenes, "FSM_NNC_class"] <- "Genes expressing\nboth NNC and\n FSM isoforms"
    isoPerGene[isoPerGene$associatedGene %in% NNC_just_genes, "FSM_NNC_class"] <- "Genes expressing\n only NNC isoforms"
    isoPerGene[isoPerGene$associatedGene %in% FSM_just_genes, "FSM_NNC_class"] <- "Genes expressing\n only FSM isoforms"
    sqantiData[sqantiData$associated_gene %in% FSMandNNCgenes, "class"] <- "Genes expressing\nboth NNC and\n FSM isoforms"
    sqantiData[sqantiData$associated_gene %in% NNC_just_genes, "class"] <- "Genes expressing\n only NNC isoforms"
    sqantiData[sqantiData$associated_gene %in% FSM_just_genes, "class"] <- "Genes expressing\n only FSM isoforms"
    
    isoPerGene$FSM_NNC_class = factor(isoPerGene$FSM_NNC_class, levels=c("Genes expressing\nboth NNC and\n FSM isoforms","Genes expressing\n only NNC isoforms","Genes expressing\n only FSM isoforms"),
                                labels=c("Genes expressing\nboth NNC and\n FSM isoforms","Genes expressing\n only NNC isoforms","Genes expressing\n only FSM isoforms"), order=T)
    
    p13 <- ggplot(data=isoPerGene[!is.na(isoPerGene$FSM_NNC_class),], aes(x=FSM_NNC_class, y=log2(geneExp+1), fill=FSM_NNC_class)) +
      geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
      ylab("log2( # Short reads per gene + 1)") +
      theme(axis.title.x=element_blank()) +
      #theme(plot.margin = unit(c(1.5,1,0.5,1), "cm")) +
      scale_fill_manual(values = c("grey38",myPalette[[4]],myPalette[[1]])) +
      guides(fill=FALSE) +
      mytheme +
      theme(axis.title.x=element_blank()) +
      ggtitle("Gene expression level in NNC/FSM containing genes\n\n" ) +
      scale_x_discrete(breaks=c("Genes expressing\nboth NNC and\n FSM isoforms",
                                "Genes expressing\n only FSM isoforms",
                                "Genes expressing\n only NNC isoforms"),
                       labels=c("NNC/FSM genes",
                                "FSM genes",
                                "NNC genes"), drop=FALSE) 
    
    
    p13.c <- ggplot(data=sqantiData[!is.na(sqantiData$class),], aes(x=class, y=log2(iso_exp+1), fill=structural_category)) +
      geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
      ylab("log2( # Short reads per transcript + 1)") +
      theme(axis.title.x=element_blank()) +
      #theme(plot.margin = unit(c(1.5,1,0.5,1), "cm")) +
      scale_fill_manual(values = myPalette) +
      guides(fill=FALSE) +
      mytheme +
      theme(axis.title.x=element_blank()) +
      ggtitle("Transcript expression level in NNC/FSM containing genes\n\n" ) +
      scale_x_discrete(breaks=c("Genes expressing\nboth NNC and\n FSM isoforms",
                                "Genes expressing\n only FSM isoforms",
                                "Genes expressing\n only NNC isoforms"),
                       labels=c("NNC/FSM genes",
                                "FSM genes",
                                "NNC genes"), drop=F) 
    
  }
}



# PLOT 23: Junction categories

if (nrow(junctionsData) > 0){
    
    junctionsData$junctionLabel = with(junctionsData, paste(chrom, strand,genomic_start_coord, genomic_end_coord, sep="_"))
  
    junctionsData$canonical_known = with(junctionsData, paste(junction_category,canonical,"SJ", sep="_"))
    junctionsData$canonical_known=as.factor(junctionsData$canonical_known)
    junctionsData$canonical_known = factor(junctionsData$canonical_known, levels=c("known_canonical_SJ", "known_non_canonical_SJ", "novel_canonical_SJ", "novel_non_canonical_SJ"),
                                           labels=c("Known\ncanonical ", "Known\nNon-canonical ", "Novel\ncanonical ", "Novel\nNon-canonical "), order=T) 
    junctionsData$structural_category = sqantiData[junctionsData$isoform,"structural_category"]
    junctionsData$TSSrange =cut(junctionsData$transcript_coord, breaks = c(0,40,80,120,160,200,10000000), labels = c("0-40", "41-80", "81-120", "121-160", "161-200",">200"))
    
    p23 <- ggplot(data=junctionsData, aes(x=structural_category)) +
      geom_bar(position="fill", aes(y = (..count..)/sum(..count..), fill=canonical_known), color="black",  size=0.3, width = 0.7) +
      scale_y_continuous(labels = percent, expand = c(0,0)) +
      scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=FALSE) +
      ylab("% Splice junctions") +
      mytheme +
      guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.3))+
      theme(legend.position="bottom", legend.title=element_blank())  +
      theme(axis.text.x = element_text(angle = 45)) +
      theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12))+
      theme(axis.title.x=element_blank()) +
      #theme(plot.margin = unit(c(1.5,1,0,1), "cm")) +
      ggtitle("Distribution of SJ type among structural classification\n\n\n")
    
    sqantiData$AllCanonical = factor(sqantiData$all_canonical, 
                                     levels = c("canonical","non_canonical"), 
                                     ordered=TRUE)
    
    p23.b <- ggplot(data=sqantiData[which(sqantiData$exons>1),], aes(x=structural_category)) +
      geom_bar(position="fill", aes(y = (..count..)/sum(..count..), fill=AllCanonical), color="black", size=0.3, width = 0.7) +
      scale_y_continuous(labels = percent, expand = c(0,0)) +
      scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=FALSE) +
      ylab("% Transcripts ") +  
      mytheme +
      guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.3))+
      theme(legend.position="bottom", legend.title=element_blank())  +
      theme(axis.text.x = element_text(angle = 45)) + 
      theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12))+
      theme(axis.title.x=element_blank()) + 
      #theme(plot.margin = unit(c(1.5,1,0,1), "cm")) +
      labs(   title="Distribution of transcripts by splice junction category\n\n\n",
              subtitle="Non canonical transcripts are those with at least one non-canonical junction\n\n") 
    
    #alpha
    # p23.c <- ggplot(data=sqantiData, aes(x=structural_category, alpha=all_canonical, fill=structural_category)) +
    #   geom_bar(position="fill", aes(y = (..count..)/sum(..count..)), color="black",  size=0.3,width = 0.7) +
    #   scale_fill_manual(values = myPalette, guide='none' )  +
    #   scale_y_continuous(labels = percent, expand = c(0,0)) +
    #   scale_alpha_manual(values=c(1,0.3)) +
    #   ylab("% Transcripts ") +  
    #   mytheme +
    #   theme(legend.position="bottom", axis.title.x = element_blank()) +
    #   theme(axis.text.x = element_text(angle = 45)) + 
    #   theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12))+
    #   theme(axis.title.x=element_blank()) + 
    #   #theme(plot.margin = unit(c(1.5,1,0,1), "cm")) +
    #   labs(   title="Distribution of transcripts by splice junction category\n\n\n",
    #           subtitle="Non canonical transcripts are those with at least one non-canonical junction\n\n") 
    # 
}


# PLOT 24-26: Junction distance to TSS

if (nrow(junctionsData) > 0){
  
  p24 <- ggplot(data=junctionsData, aes(x=transcript_coord, fill = canonical_known)) +
    geom_density(alpha=0.7,  size=0.3) +
    scale_y_continuous(expand = c(0,0))  +
    scale_x_continuous(expand = c(0,0)) +
    scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=FALSE) +
    #geom_vline(xintercept = 200, linetype = "longdash", color="red") +
    mytheme +
    guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.7))+
    theme(legend.position="bottom" ) +
    labs(list(fill = "Junction Type", x= "Distance to TSS (bp)", 
              title = "Distribution of splice junctions distance to TSS\n\n") )

  p25 <- ggplot(data=junctionsData, aes(x=TSSrange, fill=canonical_known)) +
    geom_bar(aes(y = (..count..)/sum(..count..)), color="black", size=0.3, width=0.7, position="fill") +
    scale_y_continuous(labels = percent, expand = c(0,0)) +
    scale_x_discrete(drop=FALSE) +
    scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=F) +
    mytheme +
    theme(legend.position="bottom") +
    guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.7))+
    labs(list(fill = "Junction Type", x= "TSS distance range (bp)", y="% Junctions",
              title = "Splice junction distance to TSS across junction type\n\n\n") ) +
    theme(axis.title.x  = element_text(margin=margin(10,0,0,0), size=12))

  p26 <- ggplot(data=junctionsData, aes(x=canonical_known, fill=TSSrange)) +
    geom_bar(aes(y = (..count..)/sum(..count..)), color="black", size=0.3, width=0.7, position="fill") +
    scale_y_continuous(labels = percent, expand = c(0,0)) +
    scale_fill_manual(values = myPalette, drop=F) +
    scale_x_discrete(drop=FALSE) +
    ylab("% Junctions") +
    mytheme +
    theme(legend.position="bottom") +
    guides(fill = guide_legend(title = "TSS distance range (bp)")) +
    ggtitle( "Splice junction distance to TSS across junction type\n\n\n") +
    theme(axis.title.x=element_blank()) 
  
}


# PLOT 29: RT-switching

if (nrow(junctionsData) > 0){
  
  if (nrow(junctionsData[junctionsData$RTS_junction=="TRUE",])!=0){
    
    a = data.frame(table(junctionsData$canonical_known))
    b = data.frame(table(junctionsData[which(junctionsData$RTS_junction==TRUE),"canonical_known"]))
    
    a_b_Rts = merge(a,b, by="Var1")
    a_b_Rts$perc = a_b_Rts$Freq.y/a_b_Rts$Freq.x *100
    a_b_Rts[is.na(a_b_Rts$perc),"perc"] <- 0
    
    max_height_rts = max(a_b_Rts$perc)
    mm_rts = roundUp(max_height_rts)/10
    max_height_rounded_rts = RoundUp(max_height_rts, mm_rts)
    
    p29 <- ggplot(data=a_b_Rts, aes(x=Var1, y=perc, fill=Var1)) +
      geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
      scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=F) +
      ylab("% RT-switching junctions") +
      ggtitle( "RT-switching by splice junction category \n\n" ) +
      mytheme +
      guides(fill=FALSE) +
      scale_y_continuous(expand = c(0,0), limits = c(0,max_height_rounded_rts)) +
      theme(axis.text.x = element_text(size=11))
    
    
    uniqJuncRTS =unique(junctionsData[,c("junctionLabel","canonical_known", "RTS_junction")])
    
    c = data.frame(table(uniqJuncRTS$canonical_known))
    d = data.frame(table(uniqJuncRTS[which(uniqJuncRTS$RTS_junction==TRUE),"canonical_known"]))
    
    c_d_Rts = merge(c,d, by="Var1")
    c_d_Rts$perc = c_d_Rts$Freq.y/c_d_Rts$Freq.x *100
    c_d_Rts[is.na(c_d_Rts$perc),"perc"] <- 0
    
    max_height_rts = max(c_d_Rts$perc)
    mm_rts = roundUp(max_height_rts)/10
    max_height_rounded_rts = RoundUp(max_height_rts, mm_rts)
    
    p29.a <- ggplot(data=c_d_Rts, aes(x=Var1, y=perc, fill=Var1)) +
      geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
      scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=FALSE) +
      ylab("% RT-switching junctions") +
      ggtitle( "RT-switching by splice junction category (unique junctions) \n\n" ) +
      mytheme +
      guides(fill=FALSE) +
      scale_y_continuous(expand = c(0,0), limits = c(0,max_height_rounded_rts)) +
      theme(axis.text.x = element_text(size=11))
    
  }
}

# PLOT pn4-5: Splice Junction Coverage (if coverage provided)

if (nrow(junctionsData) > 0){
    
  if (!all(is.na(junctionsData$total_coverage))){
    
    uniqJuncCov =unique(junctionsData[,c("junctionLabel","canonical_known", "total_coverage")])
    
    e = data.frame(table(uniqJuncCov$canonical_known))
    f = data.frame(table(uniqJuncCov[which(uniqJuncCov$total_coverage==0),"canonical_known"]))
    
    e_f_Rts = merge(e,f, by="Var1")
    e_f_Rts$perc = e_f_Rts$Freq.y/e_f_Rts$Freq.x *100
    e_f_Rts[is.na(e_f_Rts$perc),"perc"] <- 0
    
    max_height_rts = max(e_f_Rts$perc)
    mm_rts = roundUp(max_height_rts)/10
    max_height_rounded_rts = RoundUp(max_height_rts, mm_rts)
    
    pn4 <-ggplot(data=e_f_Rts, aes(x=Var1,fill=Var1, y=Freq.y)) +
      geom_bar(stat="identity", position = position_dodge(), color="black", size=0.3, width=0.7) +
      scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=FALSE) +
      scale_y_continuous( expand = c(0,0)) +
      ylab("# Splice Junctions without coverage") + 
      xlab("Junction Type") +
      mytheme +
      guides(fill=FALSE) +
      ggtitle( "Splice junctions without short-read coverage (unique junctions) \n\n\n") 
    
    
    pn5 <- ggplot(data=e_f_Rts, aes(x=Var1, y=perc, fill=Var1)) +
      geom_bar(position=position_dodge(), stat="identity", color="black", size=0.3, width=0.7) +
      guides(fill=FALSE) +
      scale_y_continuous( expand = c(0,0)) +
      scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=FALSE) +
      ylab("% Splice Junctions without coverage") + 
      xlab("Junction Type") +
      mytheme +
      guides(fill=FALSE) +
      ggtitle( "Splice junctions without short-read coverage (unique junctions) \n\n\n") 
    
    
  }
}
  
  
# PLOT pn1.2: Splice Junction relative coverage (if coverage and expression provided)

if (nrow(junctionsData) > 0){
  
  if (!all(is.na(junctionsData$total_coverage)) & !all(is.na(sqantiData$iso_exp))){
    
    junctionsData$isoExp = sqantiData[junctionsData$isoform, "iso_exp"]
  
    total = aggregate(cbind(total_coverage,isoExp,transcript_coord) ~ junctionLabel, data = junctionsData, 
                      FUN = function(x) c(mn = sum(x), n = min(x) ) )
    
    total$relCov = total$total_coverage[,"n"] / total$isoExp[,"mn"]
    total$minTSS = total$transcript_coord[,"n"]
    
    uniqJunc = unique(junctionsData[,c("junctionLabel", "canonical_known", "total_coverage")])
    uniqJunc$notCov = uniqJunc$total_coverage == 0
    
    uniqueJunc_nonCov = as.data.frame(table(uniqJunc[uniqJunc$totalCoverage==0,"canonical_known"])/table(uniqJunc$canonical_known)*100)
    
    uniqJunc2 = merge(total, uniqJunc, by=1)
    uniqJunc2$TSSrange =cut(uniqJunc2$minTSS, breaks = c(0,40,80,120,160,200,10000000), labels = c("0-40", "41-80", "81-120", "121-160", "161-200",">200"))
    


    # calculate total expression associated to each unique junction
    sumExpPerJunc = tapply(junctionsData$isoExp, junctionsData$junctionLabel, sum)  
    
    junctionsData$sumIsoExp = sumExpPerJunc[junctionsData$junctionLabel]
    
    junctionsData$relCov = junctionsData$total_coverage / junctionsData$sumIsoExp
    
    max_dist = max(junctionsData$transcript_coord) +1
    
    junctionsData$TSSrange = cut(junctionsData$transcript_coord, breaks = c(0,20,40,60,80,100,120,140,160,180,200,max_dist), labels = c("0-20", "21-40","41-80","61-80", "81-100","101-120", "121-140","141-160", "161-180", "181-200", ">200"))
    
    pn1.2 <-ggplot(data=junctionsData[junctionsData$relCov<1,], aes(y=relCov,x=TSSrange,fill=canonical_known)) +
      geom_boxplot(outlier.size = 0.2, size=0.3) +
      scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=FALSE) +
      ylab("Relative coverage") + 
      xlab("# TSS distance range") +
      mytheme_bw +
      theme(legend.position="bottom", legend.title=element_blank())  +
      ggtitle( "Relative Coverage of junctions\n\n\n") +
      theme(axis.text.x = element_text(angle = 45,margin=margin(15,0,0,0), size=12)) 
    
    
  }else{    uniqJunc = unique(junctionsData[,c("junctionLabel", "canonical_known")])
}

}




# PLOT p21-22: Full-lengthness (if FSM transcripts)

if (nrow(sqantiData[sqantiData$structural_category=="FSM",])!=0){

  perfect_match = sqantiData[which(sqantiData$structural_category=="FSM"), 
                             c("isoform", "length", "exons", "associated_gene", "ref_length", "ref_exons" ,"diff_to_TSS", "diff_to_TTS")]
  
  perfect_match = perfect_match[which(perfect_match$exons!=1),]

  max = max(c(abs(perfect_match$diff_to_TTS), abs(perfect_match$diff_to_TSS))) 
  perfect_match$rangeDiffTTS = cut(-(perfect_match$diff_to_TTS), breaks = c((max+1)*(-1),seq(-200, 200, by = 20),max+1))
  perfect_match$rangeDiffTSS = cut(-(perfect_match$diff_to_TSS), breaks = c((max+1)*(-1),seq(-200, 200, by = 20),(max+1)))
  
  max_height = max(c(max(table(perfect_match$rangeDiffTSS)),max(table(perfect_match$rangeDiffTTS))))
  mm = roundUp(max_height)/10
  max_height_rounded = RoundUp(max_height, mm)
  
  p21.bbbb <- ggplot(data=perfect_match, aes(x=rangeDiffTTS)) +
    geom_bar(fill=myPalette[4], color="black", size=0.3) +
    scale_y_continuous(expand = c(0,0), limits = c(0,max_height_rounded))+
    mytheme +
    scale_x_discrete(drop=F) +
    ylab("# FSM transcripts")+
    xlab("Distance to annotated TTS (bp)")+
    labs(     title="Distance distribution from sequenced to annotated TTS \n\n",
              subtitle="Negative values indicate that the sequenced TTS is upstream annotated TTS\n\n") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  p21.bbbbb <- ggplot(data=perfect_match, aes(x=rangeDiffTTS)) +
    geom_bar(aes(y = (..count..)/sum(..count..)), fill=myPalette[4], color="black", size=0.3)+
    scale_y_continuous(labels = percent_format(), limits = c(0,1), expand = c(0,0))+
    scale_x_discrete(drop=F) +
    mytheme +
    ylab("% FSM transcripts")+
    xlab("Distance to annotated TTS (bp)")+
    labs(     title="Distance distribution from sequenced to annotated TTS \n\n",
              subtitle="Negative values indicate that the sequenced TTS is upstream annotated TTS\n\n") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  
  p22.bbbbb <- ggplot(data=perfect_match, aes(x=rangeDiffTSS)) +
    geom_bar(aes(y = (..count..)/sum(..count..)), fill=myPalette[6], color="black", size=0.3)+
    scale_y_continuous(labels = percent_format(), limits = c(0,1), expand = c(0,0))+
    scale_x_discrete(drop=F) +
    mytheme +
    ylab("% FSM transcripts")+
    xlab("Distance to annotated TSS (bp)")+
    labs(     title="Distance distribution from sequenced to annotated TSS \n\n",
              subtitle="Negative values indicate that the sequenced TSS is downstream annotated TSS\n\n") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  
  p22.bbbb <- ggplot(data=perfect_match, aes(x=rangeDiffTSS)) +
    geom_bar(fill=myPalette[6], color="black", size=0.3)+
    scale_y_continuous(expand = c(0,0), limits = c(0,max_height_rounded))+
    scale_x_discrete(drop=F) +
    mytheme +
    ylab("# FSM transcripts")+
    xlab("Distance to annotated TSS (bp)")+
    labs(     title="Distance distribution from sequenced to annotated TSS \n\n",
              subtitle="Negative values indicate that the sequenced TSS is downstream annotated TSS\n\n") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  
}
  

# PLOT p28: Attribute summary if junctions

if (nrow(junctionsData) > 0){
  
  if (!all(is.na(junctionsData$total_coverage))){
    
    sqantiData[which(sqantiData$min_cov==0), "Not coverage SJ"] <- "Not coverage SJ"
    a = sqantiData[,c("all_canonical", "Not coverage SJ", "RTS_stage", "isoform", "structural_category")]
    a$RTS_stage = as.factor(a$RTS_stage)
    b = melt(a, id.vars = c("isoform", "structural_category"), na.rm = T)
    b = b[-which(b$value=="FALSE" | b$value=="canonical"),]
    b$value <- factor(b$value)
    b$value = factor(b$value, 
                     labels = c("Non-canonical SJ","RT-switching", "Not coverage SJ"), 
                     levels = c("non_canonical", "TRUE", "Not coverage SJ"), 
                     ordered=TRUE)

    
    numberPerCategory = as.data.frame(table(sqantiData$structural_category))
    attribute.df2 = merge(numberPerCategory, as.data.frame(table(b$structural_category, b$value)), by="Var1")
    attribute.df2$perc = attribute.df2$Freq.y / attribute.df2$Freq.x * 100
    attribute.df2[is.na(attribute.df2$perc),"perc"] <- 0

    p28.a <- ggplot(data=attribute.df2[attribute.df2$Var1%in%c("FSM", "NNC", "NIC"),], aes(x=Var1, y=perc, fill= Var2)) +
      geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
      scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
      scale_fill_manual(values = myPalette[9:11]) +
      ylab("% Transcripts") + 
      xlab("") +
      mytheme +
      theme(legend.position="bottom", axis.title.x = element_blank()) +
      ggtitle( "Quality control attributes across structural categories\n\n" ) +
      guides(fill = guide_legend(title = "QC Attributes") )
    
  }else{
   
    a = sqantiData[,c("all_canonical", "RTS_stage", "isoform", "structural_category")]
    a$RTS_stage = as.factor(a$RTS_stage)
    b = melt(a, id.vars = c("isoform", "structural_category"), na.rm = T)
    numberPerCategory = as.data.frame(table(sqantiData$structural_category))
    b = b[-which(b$value=="FALSE" | b$value=="canonical"),]
    b$value <- factor(b$value)
    b$value = factor(b$value, 
                     labels = c("Non-canonical SJ","RT-switching"), 
                     levels = c("non_canonical", "TRUE"), ordered=T)
    attribute.df2 = merge(numberPerCategory, as.data.frame(table(b$structural_category, b$value)), by="Var1")
    attribute.df2$perc = attribute.df2$Freq.y / attribute.df2$Freq.x * 100
    attribute.df2[is.na(attribute.df2$perc),"perc"] <- 0

    p28.a <- ggplot(data=attribute.df2[attribute.df2$Var1%in%c("FSM", "NNC", "NIC"),], aes(x=Var1, y=perc, fill= Var2)) +
      geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
      scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
      scale_fill_manual(values = myPalette[9:11]) +
      ylab("% Transcripts") + 
      xlab("") +
      mytheme +
      theme(legend.position="bottom", axis.title.x = element_blank()) +
      ggtitle( "Quality control attributes across structural categories\n\n" ) +
      guides(fill = guide_legend(title = "QC Attributes") )    
    
    
    
  }
    }


# PLOT p30,p31,p32: percA by subcategory


p30 <- ggplot(data=sqantiData, aes(y=perc_A_downstream_TTS, x=structural_category, fill=subcategory)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) + mytheme +
  ylab("% Adenines in gDNA window downstreamTTS") +
  mytheme  + theme(axis.text.x = element_text(angle = 45)) +
  theme(legend.position="bottom", legend.title=element_blank(), legend.direction = "horizontal", legend.box = "vertical") +
  guides(fill=guide_legend(nrow=5,byrow=TRUE)) +
  theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12))+
  ggtitle("Adenine enrichment downstream TTS (intra primming)\n\n" ) +
  theme(axis.title.x=element_blank()) +
  scale_fill_manual(values=myPalette, breaks=c("3prime_fragment", "internal_fragment", "5prime_fragment",
                            "mono-exon", "multi-exon", "combination_of_known_junctions", 
                            "no_combination_of_known_junctions", "mono-exon_by_intron_retention/s",
                            "not any annotated donor/acceptor", "any annotated donor/acceptor"),
                   labels=c("3' fragment", "Internal fragment", "5' fragment",
                            "Mono-exon", "Multi-exon", "Combination of annotated junctions", 
                            "Not combination of annotated junctions", "Mono-exon by intron retention",
                            "Without annotated donors/acceptors", "At least one annotated donor/acceptor"), drop=F) 

p31 <- ggplot(data=sqantiData, aes(y=perc_A_downstream_TTS, x=structural_category, fill=exonCat)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) + mytheme +
  scale_fill_manual(breaks=c("MonoExon", "MultiExon"),
                 labels=c("Mono-exon Isoforms", "Multi-exon Isoforms"), values=myPalette) +
  ylab("% Adenines in gDNA window downstreamTTS") +
  theme(legend.position="bottom", legend.title=element_blank()) +
  theme(axis.text.x = element_text(angle = 45)) + 
  theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12))+
  ggtitle("Adenine enrichment downstream TTS (intra primming)\n\n" ) +  
  theme(axis.title.x=element_blank()) 


p32 <- ggplot(data=sqantiData, aes(y=perc_A_downstream_TTS, x=structural_category, fill=coding)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) + mytheme +
  scale_fill_manual(breaks=c("coding", "non_coding"),
                    labels=c("Coding Isoforms", "Non-coding Isoforms"), values=myPalette[3:4]) +
  ylab("% Adenines in gDNA window downstreamTTS") +
  theme(legend.position="bottom", legend.title=element_blank() ) +
  theme(axis.text.x = element_text(angle = 45)) + 
  theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12))+
  ggtitle("Adenine enrichment downstream TTS (intra primming)\n\n" ) +  
  theme(axis.title.x=element_blank()) 









###** Output plots

pdf(file=report.file, width = 6.5, height = 6.5)


#cover
grid.newpage()
cover <- textGrob("SQANTI report",
                  gp=gpar( fontface="italic", fontsize=40, col="orangered"))
grid.draw(cover)

# t1
freqCat = as.data.frame(table(sqantiData$structural_category))
freqCat = freqCat[order(freqCat$Freq,decreasing = T),]
table <- tableGrob(freqCat, rows = NULL, cols = c("category","# isoforms"))
title <- textGrob("Characterization of transcripts\n based on splice junctions", gp=gpar(fontface="italic", fontsize=17), vjust = -3.5)
gt1 <- gTree(children=gList(table, title))

# t2
freqCat = as.data.frame(table(isoPerGene$novelGene))
table2 <- tableGrob(freqCat, rows = NULL, cols = c("category","# genes"))
title2 <- textGrob("Gene classification", gp=gpar(fontface="italic", fontsize=17), vjust = -3.5)
gt2 <- gTree(children=gList(table2, title2))

# t4
if (nrow(junctionsData) > 0){
  freqCat = as.data.frame(table(uniqJunc$canonical_known))
  freqCat$Var1 = gsub(" ", "", freqCat$Var1)
  freqCat$Var1 = gsub("\n", " ", freqCat$Var1)
  table2 <- tableGrob(freqCat, rows = NULL, cols = c("category","# SJ"))
  title2 <- textGrob("SJ classification", gp=gpar(fontface="italic", fontsize=17), vjust = -5)
  gt3 <- gTree(children=gList(table2, title2))
}else{  
  title2 <- textGrob("SJ classification", gp=gpar(fontface="italic", fontsize=17), vjust = -5)
  freqCat = data.frame(Var1=c("Known canonical", "Known Non-canonical", "Novel canonical", "Novel Non-canonical"), Freq=c(0,0,0,0))
  table2 <- tableGrob(freqCat, rows = NULL, cols = c("category","# SJ"))
  gt3 <- gTree(children=gList(table2, title2))}

# t3
nGenes = nrow(isoPerGene)
nIso = nrow(sqantiData)
sn = paste("# Genes: ", nGenes, "\n", "# Isoforms: ", nIso)
s <- textGrob(sn, gp=gpar(fontface="italic", fontsize=17), vjust = 0)


# drawing t1 and t2
grid.arrange(s,gt2,gt3,gt1, layout_matrix = cbind(c(1,2,3),c(1,4,4)))


#1. general parameters
s <- textGrob("Gene characterization", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
grid.arrange(s)
p0
p7
p6
if (!all(is.na(sqantiData$iso_exp))){
  print(p10)
}
if (!all(is.na(sqantiData$FL))){
  print(p11)
}
print(p31)
print(p32)

#2. general parameters by structual categories

s <- textGrob("Structrual Isoform characterization\nbased on splice junctions", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
grid.arrange(s)
print(p1)
print(p4)
print(p5)
if (!all(is.na(sqantiData$iso_exp))){
  print(p8)
}
if (!all(is.na(sqantiData$FL))){
  print(p9)
}
print(p30)

if (nrow(sqantiData[sqantiData$structural_category%in%c("FSM","ISM"),])!=0){
  print(p2)
  print(p3)
}
if (!all(is.na(sqantiData$gene_exp))){
  if (nrow(sqantiData[sqantiData$structural_category=="NNC",])!=0){
    print(p12)
  }
}
if (!all(is.na(sqantiData$gene_exp))){
  if (nrow(sqantiData[sqantiData$structural_category=="NNC",])!=0 & nrow(sqantiData[sqantiData$structural_category=="FSM",])!=0 ){
    print(p13)
    print(p13.c)
  }
}


#3. splice junction

if (nrow(junctionsData) > 0){
 
  s <- textGrob("Splice junction characterization", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
  grid.arrange(s)
  print(p23.b)
  print(p23)
  print(p24)
  print(p25)
  print(p26)
  
  if (!all(is.na(junctionsData$total_coverage)) & !all(is.na(sqantiData$iso_exp))){
    print(pn1.2)
  }
  
  if (!all(is.na(junctionsData$total_coverage))){
    print(pn4)
    print(pn5)
  }
  
  if (nrow(junctionsData[junctionsData$RTS_junction=="TRUE",])!=0){
    print(p29)
    print(p29.a)
  }
} else{
  s <- textGrob("Splice junction characterization can be plotted:  Not junctions found in input transcriptome", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
  grid.arrange(s)
}


#4. full-lengthness

if (nrow(sqantiData[sqantiData$structural_category=="FSM",])!=0){
  s <- textGrob("Full-lengthness characterization of isoforms", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
  grid.arrange(s)
  print(p21.bbbb)
  print(p21.bbbbb)
  print(p22.bbbb)
  print(p22.bbbbb)
  
}

if (nrow(junctionsData) > 0){
  s <- textGrob("Quality control attributes", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
  grid.arrange(s)
  print(p28.a)
} else{  s <- textGrob("Quality control attributes can be computed: Not junctions found in input transcriptome", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
grid.arrange(s)}



dev.off()


print("SQANTI report successfully generated!")

















