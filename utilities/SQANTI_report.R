###############################
##### SQANTI report generation
###############################

### Author: Lorena de la Fuente
### Date: 14/10/2016

#********************** Taking argument from python script

args <- commandArgs(trailingOnly = TRUE)
class.file = args[1]
junc.file = args[2]
report.file = paste(strsplit(class.file, "_classification.txt")[[1]][1], "Report.pdf", sep="_")


#********************** Packages (installed if not found)

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

sorted <- sqantiData[order(sqantiData$isoExp, decreasing = T),]
FSMhighestExpIsoPerGene <- sorted[(!duplicated(sorted$associatedGene) & sorted$structuralCategory=="full-splice_match"),"isoform"]
sqantiData[which(sqantiData$isoform%in%FSMhighestExpIsoPerGene),"RTS_stage"] <- FALSE
write.table(sqantiData, file=class.file, row.names=FALSE, quote=F, sep="\t")

xaxislabelsF1 = c("FSM", "ISM", "NIC", "NNC", "Genic\nGenomic",  "Antisense", "Fusion","Intergenic", "Genic\nIntron")
sqantiData$structuralCategory = factor(sqantiData$structuralCategory, 
                                       labels = xaxislabelsF1, 
                                       levels = c("full-splice_match","incomplete-splice_match","novel_in_catalog","novel_not_in_catalog", "genic","antisense","fusion","intergenic","genic_intron"), 
                                       ordered=TRUE)




########### Junction information

junctionsData = read.table(file=junc.file, header=T, as.is=T, sep="\t")
junctionsData[which(junctionsData$isoform%in%FSMhighestExpIsoPerGene),"RTS_junction"] <- FALSE
write.table(junctionsData, file=junc.file, row.names=FALSE, quote=F, sep="\t")


########### Handling information 

sqantiData[grep("novelGene",  sqantiData$associatedGene), "novelGene"] <- "Novel Genes"
sqantiData[is.na(sqantiData$novelGene), "novelGene"] <- "Annotated Genes"
sqantiData[which(sqantiData$exons>1), "exonCat"] <- "MultiExon"
sqantiData[-which(sqantiData$exons>1), "exonCat"] <- "MonoExon"
#sorted <- sqantiData[order(sqantiData$isoExp, decreasing = T),]
#sorted$highestExpIso <- !duplicated(sorted$associatedGene)
#sqantiData = sorted
#sqantiData[which(sqantiData$highestExpIso==TRUE & sqantiData$structuralCategory=="FSM"),"RTS_stage"] <- "FALSE"

junctionsData$canonical_known = with(junctionsData, paste(junctionCategory,canonical,"SJ", sep="_"))
junctionsData$canonical_known=as.factor(junctionsData$canonical_known)
junctionsData$canonical_known = factor(junctionsData$canonical_known, levels=c("known_canonical_SJ", "known_non_canonical_SJ", "novel_canonical_SJ", "novel_non_canonical_SJ"),
                                       labels=c("Known\ncanonical ", "Known\nNon-canonical ", "Novel\ncanonical ", "Novel\nNon-canonical "), order=T) 
junctionsData$structuralCategory = sqantiData[junctionsData$isoform,"structuralCategory"]
junctionsData$TSSrange =cut(junctionsData$transcriptCoord, breaks = c(0,40,80,120,160,200,10000000), labels = c("0-40", "41-80", "81-120", "121-160", "161-200",">200"))



# Gene level information
isoPerGene = aggregate(sqantiData$isoform, by = list("associatedGene" = sqantiData$associatedGene, "novelGene" = sqantiData$novelGene, "FSM_class" = sqantiData$FSM_class, "geneExp"=sqantiData$geneExp), length)
colnames(isoPerGene)[ncol(isoPerGene)] <- "nIso"

isoPerGene$range =cut(isoPerGene$nIso, breaks = c(0,1,3,5,2500), labels = c("1", "2-3", "4-5", ">=6"))

FL_gene = aggregate(as.integer(sqantiData$FL), by = list("associatedGene" = sqantiData$associatedGene), sum)
colnames(FL_gene)[ncol(FL_gene)] <- "FL_gene"
gene = merge(isoPerGene, FL_gene, by="associatedGene")
gene$FSM_class2 = factor(gene$FSM_class, 
                         levels = c("A", "B", "C"), 
                         labels = c("MonoIsoform Gene", "MultiIsoform Genes\nwithout expression\nof a FSM", "MultiIsoform Genes\nexpressing at least\none FSM"), 
                         ordered=TRUE)


# Genes expression only FSM, only NNC and both
FSM_just_genes = unique(sqantiData[sqantiData$FSM_class=="A" & sqantiData$structuralCategory=="FSM","associatedGene"]) 
NNC_just_genes = unique(sqantiData[sqantiData$FSM_class=="A" & sqantiData$structuralCategory=="NNC","associatedGene"]) 
FSMandNNCgenes = unique(sqantiData[sqantiData$FSM_class=="C" & sqantiData$structuralCategory=="NNC","associatedGene"]) 
gene[gene$associatedGene %in% FSMandNNCgenes, "FSM_NNC_class"] <- "Genes expressing\nboth NNC and\n FSM isoforms"
gene[gene$associatedGene %in% NNC_just_genes, "FSM_NNC_class"] <- "Genes expressing\n only NNC isoforms"
gene[gene$associatedGene %in% FSM_just_genes, "FSM_NNC_class"] <- "Genes expressing\n only FSM isoforms"
sqantiData[sqantiData$associatedGene %in% FSMandNNCgenes, "class"] <- "Genes expressing\nboth NNC and\n FSM isoforms"
sqantiData[sqantiData$associatedGene %in% NNC_just_genes, "class"] <- "Genes expressing\n only NNC isoforms"
sqantiData[sqantiData$associatedGene %in% FSM_just_genes, "class"] <- "Genes expressing\n only FSM isoforms"

gene$FSM_NNC_class = factor(gene$FSM_NNC_class, levels=c("Genes expressing\nboth NNC and\n FSM isoforms","Genes expressing\n only NNC isoforms","Genes expressing\n only FSM isoforms"),
                                       labels=c("Genes expressing\nboth NNC and\n FSM isoforms","Genes expressing\n only NNC isoforms","Genes expressing\n only FSM isoforms"), order=T)

# NNC expression genes vs not NNC expression genes

NNC_genes = unique(sqantiData[sqantiData$structuralCategory=="NNC","associatedGene"]) #1680
notNNC_genes = unique(sqantiData[-which(sqantiData$associatedGene%in%NNC_genes),"associatedGene"]) #6014
gene[gene$associatedGene %in% notNNC_genes, "NNC_class"] <- "Genes not expressing\n NNC isoforms"
gene[gene$associatedGene %in% NNC_genes, "NNC_class"] <- "Genes expressing\n NNC isoforms"

gene$NNC_class = factor(gene$NNC_class, levels=c("Genes expressing\n NNC isoforms","Genes not expressing\n NNC isoforms"),
                            labels=c("Genes expressing\n NNC isoforms","Genes not expressing\n NNC isoforms"), order=T)


#attribute by structural classification

sqantiData[which(sqantiData$MinCov==0), "Not coverage SJ"] <- "Not coverage SJ"
a = sqantiData[,c("AllCanonical", "Not coverage SJ", "RTS_stage", "isoform", "structuralCategory")]
a$RTS_stage = as.factor(a$RTS_stage)
a$`Not coverage SJ` = as.factor(a$`Not coverage SJ`)
b = melt(a, id.vars = c("isoform", "structuralCategory"), na.rm = T)
b = b[-which(b$value=="FALSE" | b$value=="canonical"),]
b$value <- factor(b$value)
b$value = factor(b$value, 
                 labels = c("Non-canonical SJ","RT-switching", "Not coverage SJ"), 
                 levels = c("non_canonical", "TRUE", "Not coverage SJ"), 
                 ordered=TRUE)
attribute.df = table(b$structuralCategory,b$value)
attribute.df2 = merge(as.data.frame(attribute.df), as.data.frame(table(sqantiData$structuralCategory)), by = 1)
attribute.df2$perc = attribute.df2$Freq.x / attribute.df2$Freq.y * 100



# attribute by class

sqantiData$category = with(sqantiData, paste(class,structuralCategory, sep=":"))
#sqantiData$category = factor(sqantiData$category, levels=c("Genes expressing\n NNC isoforms:NNC","Genes not expressing\n NNC isoforms"),
 #                       labels=c("Genes expressing\n NNC isoforms","Genes not expressing\n NNC isoforms"), order=T)


a.2 = sqantiData[!is.na(sqantiData$class),c("AllCanonical", "Not coverage SJ", "RTS_stage", "isoform", "category")]
a.2 = a.2[,c("AllCanonical", "Not coverage SJ", "RTS_stage", "isoform", "category")]
a.2$RTS_stage = as.factor(a.2$RTS_stage)
a.2$`Not coverage SJ` = as.factor(a.2$`Not coverage SJ`)
b.2 = melt(a.2, id.vars = c("isoform", "category"), na.rm = T)
b.2 = b.2[-which(b.2$value=="FALSE" | b.2$value=="canonical"),]
b.2$value <- factor(b.2$value)
b.2$value = factor(b.2$value, 
                   labels = c("Non-canonical SJ","RT-switching", "Not coverage SJ"), 
                   levels = c("non_canonical", "TRUE", "Not coverage SJ"), 
                   ordered=TRUE)
attribute.df.2 = table(b.2$category,b.2$value)
attribute.df.2.2 = merge(as.data.frame(attribute.df.2), as.data.frame(table(sqantiData$category)), by = 1)
attribute.df.2.2$perc = attribute.df.2.2$Freq.x / attribute.df.2.2$Freq.y * 100





# Coverage information

perfect_match = sqantiData[which(sqantiData$structuralCategory=="FSM"), 
                           c("isoform", "length", "exons", "associatedGene", "refLength", "refExons" ,"diffToTSS", "diffToTTS")]

# NEW: removing full-splice match hits that are monoexons

perfect_match = perfect_match[which(perfect_match$exons!=1),]

#5' end
perfect_match[which(perfect_match$diffToTSS<0),"end5_state"] <- "longer"
perfect_match[which(perfect_match$diffToTSS>0),"end5_state"] <- "shorter"
perfect_match[which(perfect_match$diffToTSS==0),"end5_state"] <- "equal"

perfect_match$end5_perc_coverage = (perfect_match$refLength-perfect_match$diffToTSS)/perfect_match$refLength * 100 
perfect_match$end5_perc_coverage_MAX100 = perfect_match$end5_perc_coverage
perfect_match[which(perfect_match$end5_perc_coverage>100),"end5_perc_coverage_MAX100"] <- 100
perfect_match$end5_perc_coverage_NA = perfect_match$end5_perc_coverage
perfect_match[which(perfect_match$end5_perc_coverage>100),"end5_perc_coverage_NA"] <- NA

perfect_match$range5 =cut(perfect_match$end5_perc_coverage_NA, breaks = seq(0,100,by = 5))

#3' end
perfect_match[which(perfect_match$diffToTTS<0),"end3_state"] <- "longer"
perfect_match[which(perfect_match$diffToTTS>0),"end3_state"] <- "shorter"
perfect_match[which(perfect_match$diffToTTS==0),"end3_state"] <- "equal"

perfect_match$end3_perc_coverage = (perfect_match$refLength-perfect_match$diffToTTS)/perfect_match$refLength* 100 
perfect_match$end3_perc_coverage_MAX100 = perfect_match$end3_perc_coverage
perfect_match$end3_perc_coverage_NA = perfect_match$end3_perc_coverage
perfect_match[which(perfect_match$end3_perc_coverage>(100)),"end3_perc_coverage_NA"] <- NA
perfect_match$end3_perc_coverage_MAX100 = perfect_match$end3_perc_coverage
perfect_match[which(perfect_match$end3_perc_coverage>(100)),"end3_perc_coverage_MAX100"] <- 100


q_5 = quantile(perfect_match[perfect_match$diffToTSS>=0,"diffToTSS"], probs = seq(0,1,0.1))
q_3 = quantile(perfect_match[perfect_match$diffToTTS>=0,"diffToTTS"], probs = seq(0,1,0.1))
to = data.frame(diff = c(perfect_match$diffToTSS, perfect_match$diffToTTS), type = rep(c("Nt to TSS", "Nt to TTS"), each=nrow(perfect_match)))
to = to[to$diff>=0,]
quantile95 = aggregate(to[1], to[c(2)], function(x) quantile(x,0.95))


#** Global plot parameters

#myPalette = c("#6BAED6","#FC8D59","#78C679","coral2","#969696","#66C2A4", "goldenrod1", "darksalmon", "#7BCCC4" ,"#FE9929", "#41B6C4", "#67A9CF")

myPalette = c("#6BAED6","#FC8D59","#78C679","coral2","#969696","#66C2A4", "goldenrod1", "darksalmon", "#41B6C4","tomato3", "#FE9929")

#myPalette = c("#66C2A4","#FC8D59","#6BAED6","coral2","goldenrod1","#7BCCC4","#969696", "darksalmon","#78C679" , "#FE9929", "#41B6C4", "#67A9CF")


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




#*** Generating plots


# PLOT 1

legendLabelF1 = tools::toTitleCase(gsub("_", " ", levels(as.factor(sqantiData$coding))))

p0 <- ggplot(isoPerGene, aes(x=range, fill=range)) +
  geom_bar(aes(stat = "count", y= (..count..)/sum(..count..)), color="black", size=0.3, width=0.7) +
  guides(fill=FALSE) +
  scale_y_continuous(labels = percent, expand = c(0,0)) +
  scale_fill_manual(values = myPalette[c(2:5)]) +
  labs(x ="# Isoforms per Gene", title="Distribution of isoforms per gene\n\n\n", y = "% Genes") +
  mytheme

p0.2 <- ggplot(isoPerGene, aes(x=range, fill=range)) +
  geom_bar(stat = "count",  color="black", size=0.3, width=0.7) +
  guides(fill=FALSE) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = myPalette[c(2:5)]) +
  labs(x ="# Isoforms per Gene", title="Distribution of isoforms per gene\n\n\n", y = "# Genes") +
  mytheme



p1 <- ggplot(data=sqantiData, aes(x=structuralCategory)) +
  geom_bar(aes(y = (..count..)/sum(..count..), alpha=coding, fill=structuralCategory), color="black", size=0.3, width=0.7) +
  scale_y_continuous(labels = percent, expand = c(0,0)) +
  geom_text(aes(y = ((..count..)/sum(..count..)), label = scales::percent((..count..)/sum(..count..))), stat = "count", vjust = -0.25)  +
  scale_alpha_manual(values=c(1,0.3), 
                     name = "Coding prediction", 
                     labels = legendLabelF1)+
  scale_fill_manual(values = myPalette, guide='none') + 
  xlab("") + 
  ylab("% Transcripts") +
  mytheme + 
  geom_blank(aes(y=1.05*((..count..)/sum(..count..))), stat = "count") +
  theme(axis.text.x = element_text(angle = 45)) +
  ggtitle("Isoform distribution across structural categories\n\n" ) +
  theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12)) +
  theme(legend.justification=c(1,1), legend.position=c(1,1)) 





# PLOT 2

p2 <- ggplot(data=sqantiData[sqantiData$structuralCategory%in%c("FSM","ISM"),], aes(x=structuralCategory, y=refLength, fill=structuralCategory)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) + mytheme_bw +
  scale_fill_manual(values = myPalette) +
  guides(fill=FALSE) +
  xlab("Structural Classification") +  
  ylab("Matched Reference Length (bp)") +
  ggtitle("Length distribution of matched reference transcripts\n\n" ) 





# PLOT 3

p3 <- ggplot(data=sqantiData[sqantiData$structuralCategory%in%c("FSM","ISM"),], aes(x=structuralCategory, y=refExons, fill=structuralCategory)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +   
  xlab("Structural Classification") +  
  ylab("Matched Reference exon number") +
  scale_fill_manual(values = myPalette) +
  guides(fill=FALSE) +
  mytheme_bw +
  ggtitle("Exon Number distribution of matched reference transcripts\n\n" )



# PLOT 4 

p4 <- ggplot(data=sqantiData, aes(x=structuralCategory, y=length, fill=structuralCategory)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  ylab("Transcript Length (bp)") +
  scale_fill_manual(values = myPalette) +
  guides(fill=FALSE) +
  mytheme  + theme(axis.text.x = element_text(angle = 45)) + 
  theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12))+
  ggtitle("Transcript length distribution by structural classification\n\n" ) +  
  theme(axis.title.x=element_blank()) 
#theme(plot.margin = unit(c(1.5,1,0,1), "cm")) 


# PLOT 5

p5 <- ggplot(data=sqantiData, aes(x=structuralCategory, y=exons, fill=structuralCategory)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  ylab("Number of exons per transcript") +
  scale_fill_manual(values = myPalette) +
  guides(fill=FALSE) +
  mytheme  + theme(axis.text.x = element_text(angle = 45)) + 
  theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12))+
  ggtitle("Exon number distribution by structural classification\n\n" ) +
  theme(axis.title.x=element_blank()) 
#theme(plot.margin = unit(c(1.5,1,0,1), "cm")) 


# PLOT 6

p6 <- ggplot(data=sqantiData, aes(x=novelGene)) +
  geom_bar(position="fill",aes(y = (..count..)/sum(..count..), fill=exonCat), color="black", size=0.3, width=0.5) +
  #geom_text(aes(y = (..count..), label=(..count..)), stat = "count", vjust = -0.25)  +
  scale_y_continuous(labels = percent, expand = c(0,0)) +
  scale_fill_manual(name = "Transcript type", 
                    values = myPalette[c(2:5)]) +
  ylab("% Transcripts ") +  
  mytheme +
  theme(axis.title.x=element_blank()) + 
  theme(legend.position="bottom") +   
  ggtitle("Distribution of mono/multi exon transcripts\n\n" ) 



# PLOT 7

p7 <- ggplot(data=gene, aes(x=novelGene)) +
  geom_bar(position="fill", aes(y = (..count..)/sum(..count..), fill=range), color="black", size=0.3, width=0.5) +
  scale_y_continuous(labels = percent, expand = c(0,0)) +
  scale_fill_manual(name = "# Isoforms Per Gene",
                    values = myPalette[c(2:5)]) +
  ylab("% Genes ") +  
  xlab("Gene Type") +
  mytheme +
  theme(axis.title.x=element_blank()) + 
  theme(legend.position="bottom") +   
  guides(fill = guide_legend(keywidth = 0.9, keyheight = 0.9)) +
  ggtitle("Distribution of number of isoforms\n\n\n\n" ) 
#theme(plot.margin = unit(c(1.5,1,0,1), "cm")) 




# PLOT 8

p8 <- ggplot(data=sqantiData, aes(x=structuralCategory, y=log2(isoExp+1), fill=structuralCategory)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  ylab("log2(transcript expression +1)") +
  scale_fill_manual(values = myPalette) +
  guides(fill=FALSE) +
  mytheme  + theme(axis.text.x = element_text(angle = 45)) + 
  theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12))+
  theme(axis.title.x=element_blank()) + 
  #theme(plot.margin = unit(c(1.5,1,0,1), "cm")) +
  ggtitle("Transcript expression by structural category\n\n" )


# PLOT 9 

if (!all(is.na(sqantiData$FL))){
  
  p9 <- ggplot(data=sqantiData, aes(x=structuralCategory, y=log2(FL+1), fill=structuralCategory)) +
    geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
    ylab("log2( # FL reads + 1)") +
    scale_fill_manual(values = myPalette) +
    guides(fill=FALSE) +
    mytheme +
    mytheme  + theme(axis.text.x = element_text(angle = 45)) + 
    theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12))+
    theme(axis.title.x=element_blank()) + 
    #theme(plot.margin = unit(c(1.5,1,0,1), "cm")) +
    ggtitle("Number of FL reads per transcript by structural category\n\n" )
  
}

# PLOT 10 

p10 <- ggplot(data=gene, aes(x=novelGene, y=log2(geneExp+1), fill=novelGene)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  xlab("Structural Classification") +  
  ylab("log2( # Short reads + 1)") +
  scale_fill_manual(values = myPalette[c(3:4)]) +
  guides(fill=FALSE) +
  mytheme_bw +
  theme(axis.title.x=element_blank()) + 
  #theme(plot.margin = unit(c(1.5,1,0.5,1), "cm")) +
  ggtitle("Gene expression by type of gene annotation\n\n" )


# PLOT 11

if (!all(is.na(sqantiData$FL))){
  
  p11 <- ggplot(data=gene, aes(x=novelGene, y=log2(FL_gene+1), fill=novelGene)) +
    geom_boxplot(color="black", size=0.3,outlier.size = 0.2) +
    ylab("log2( # FL reads + 1)") +
    scale_fill_manual(values = myPalette[c(3:4)]) +
    guides(fill=FALSE) +
    mytheme_bw +
    theme(axis.title.x=element_blank()) + 
    #theme(plot.margin = unit(c(1.5,1,0.5,1), "cm")) +
    ggtitle("Number of FL reads per Gene by type of gene annotation\n\n" )
  
}

# PLOT 12

p12 <- ggplot(data=gene[!is.na(gene$NNC_class),], aes(x=NNC_class, y=log2(geneExp+1), fill=NNC_class)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  xlab("") +  
  ylab("log2(# Short reads per gene + 1)") +
  scale_fill_manual(values = c(myPalette[4],"grey38")) +
  guides(fill=FALSE) +
  mytheme_bw +
  theme(axis.title.x=element_blank()) + 
  ggtitle("Gene expression levels between NNC and not NNC containing genes\n\n" ) 


# PLOT 13

p13 <- ggplot(data=gene[!is.na(gene$FSM_NNC_class),], aes(x=FSM_NNC_class, y=log2(geneExp+1), fill=FSM_NNC_class)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  ylab("log2( # Short reads per gene + 1)") +
  theme(axis.title.x=element_blank()) +
  #theme(plot.margin = unit(c(1.5,1,0.5,1), "cm")) +
  scale_fill_manual(values = c("grey38",myPalette[[4]],myPalette[[1]])) +
  guides(fill=FALSE) +
  mytheme_bw +
  theme(axis.title.x=element_blank()) +
  ggtitle("Gene Expression level in NNC/FSM containing genes\n\n" ) +
  scale_x_discrete(breaks=c("Genes expressing\nboth NNC and\n FSM isoforms",
                            "Genes expressing\n only FSM isoforms",
                            "Genes expressing\n only NNC isoforms"),
                   labels=c("NNC/FSM genes",
                            "FSM genes",
                            "NNC genes")) 


# p13.b <- ggplot(data=gene[!is.na(gene$FSM_NNC_class),], aes(x=FSM_NNC_class, y=geneExp, fill=FSM_NNC_class)) +
#   geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
#   ylab("log2( # Short reads per gene + 1)") +
#   theme(axis.title.x=element_blank()) +
#   scale_fill_manual(values = c("grey38",myPalette[[1]], myPalette[[4]])) +
#   #theme(plot.margin = unit(c(1.5,1,0.5,1), "cm")) +
#   guides(fill=FALSE) +
#   mytheme_bw +
#   theme(axis.title.x=element_blank()) +
#   ggtitle("Gene Expression level in NNC/FSM containing genes\n\n" )

p13.c <- ggplot(data=sqantiData[!is.na(sqantiData$class),], aes(x=class, y=log2(isoExp+1), fill=structuralCategory)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  ylab("log2( # Short reads per transcript + 1)") +
  theme(axis.title.x=element_blank()) +
  #theme(plot.margin = unit(c(1.5,1,0.5,1), "cm")) +
  scale_fill_manual(values = myPalette) +
  guides(fill=FALSE) +
  mytheme_bw +
  theme(axis.title.x=element_blank()) +
  ggtitle("Transcript Expression level in NNC/FSM containing genes\n\n" ) +
  scale_x_discrete(breaks=c("Genes expressing\nboth NNC and\n FSM isoforms",
                            "Genes expressing\n only FSM isoforms",
                            "Genes expressing\n only NNC isoforms"),
                   labels=c("NNC/FSM genes",
                            "FSM genes",
                            "NNC genes")) 
  
  

# p13.d <- ggplot(data=sqantiData[sqantiData$class=="Genes expressing\nboth NNC and\n FSM isoforms",], aes(y=log2(isoExp+1), x=structuralCategory, fill=structuralCategory)) +
#   geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
#   ylab("log2( # Short reads per transcript + 1)") +
#   theme(axis.title.x=element_blank()) +
#   #theme(plot.margin = unit(c(1.5,1,0.5,1), "cm")) +
#   scale_fill_manual(values = myPalette) +
#   guides(fill=FALSE) +
#   mytheme_bw +
#   theme(axis.title.x=element_blank()) +
#   ggtitle("Transcript Expression level in NNC/FSM containing genes\n\n" )


# PLOT 14

# p14 <- ggplot(data=gene, aes(x=class, y=geneExp, fill=class)) +
#   geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
#   xlab("") +
#   ylab("# Short reads per gene") +
#   scale_fill_manual(values = myPalette) +
#   guides(fill=FALSE) +
#   mytheme +
#   ggtitle("Gene Expression level for NNC and FSM containing genes\n\n" )
#

# PLOT 15

# p15 <- ggplot(data=sqantiData[!is.na(sqantiData$class),], aes(x=class, y=log2(isoExp+1), fill=structuralCategory)) +
#   geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
#   ylab("log2( # Short reads per transcript + 1)") +
#   #theme(plot.margin = unit(c(1.5,1,0.5,1), "cm")) +
#   scale_fill_manual(values = myPalette) +
#   mytheme +
#   theme(axis.title.x=element_blank()) + 
#   ggtitle("Transcript expression levels for NNC and FSM containing genes\n\n" ) +
#   guides(fill = guide_legend(title = "Structural Classification", keywidth = 0.9, keyheight = 0.9)) +
#   theme(legend.position="bottom")


# PLOT 16

# p16 <-ggplot(data=sqantiData, aes(x=class, y=log2(isoExp+1), fill=structuralCategory)) +
#   geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
#   theme(axis.title.x=element_blank()) +
#   theme(plot.margin = unit(c(1.5,1,0.5,1), "cm")) +
#   ylab("log2( # Short reads per transcripts + 1)") +
#   scale_fill_manual(values = myPalette) +
#   mytheme +
#   ggtitle("Transcript Expression distribution betweeen NNC and FSM based groups\n\n" ) +
#   guides(fill = guide_legend(title = "Structural Category", keywidth = 0.9, keyheight = 0.9)) +
#   theme(legend.position="bottom")
#   

# PLOT 17

# p17 <- ggplot() +
#   geom_density(data=perfect_match[!is.na(perfect_match$end3_perc_coverage_NA),], aes(x=end3_perc_coverage_NA), fill=myPalette[6])+
#   mytheme +
#   scale_x_continuous(limits = c(0, 100))+
#   xlab("Coverage in 3' end (%)")+
#   ylab("Density") +
#   ggtitle("Coverage of transcripts in 3' end")


# PLOT 18

# p18 <- ggplot() +
#   geom_density(data=perfect_match[!is.na(perfect_match$end5_perc_coverage_NA),], aes(x=end5_perc_coverage_NA), fill=myPalette[3])+
#   mytheme +
#   scale_x_continuous(limits = c(0, 100))+
#   xlab("Coverage in 5' end (%)")+
#   ylab("Density") +
#   ggtitle("Coverage of transcripts in 5' end")


# PLOT 19

# p19 <- ggplot(data=perfect_match[!is.na(perfect_match$end3_perc_coverage_NA) & perfect_match$end3_perc_coverage_NA>=70,], aes(x=end3_perc_coverage_NA)) +
#   geom_histogram(aes(y=..density..), fill="white", breaks=seq(70, 100, by = 2), color="black")+
#   geom_density(alpha=.6, fill=myPalette[6]) +
#   mytheme +
#   xlab("Coverage in 3' end (%)")+
#   ylab("Density") +
#   ggtitle("Coverage of transcripts in 3' end")
# 
# PLOT 20
# 
# p20 <- ggplot(data=perfect_match[!is.na(perfect_match$end5_perc_coverage_NA) & perfect_match$end5_perc_coverage_NA>=70,], aes(x=end5_perc_coverage_NA)) +
#   geom_histogram(aes(y=..density..),fill="white", breaks=seq(70, 100, by = 2), color="black") +
#   geom_density(alpha=.6, fill=myPalette[3]) +
#   mytheme +
#   xlab("Coverage in 5' end (%)")+
#   ylab("Density") +
#   ggtitle("Coverage of transcripts in 5' end")
#
# # PLOT 21
# 
# p21 <- ggplot(data=perfect_match[!is.na(perfect_match$end3_perc_coverage_NA) & perfect_match$end3_perc_coverage_NA>=70,], aes(x=end3_perc_coverage_NA)) +
#   geom_histogram(fill=myPalette[6], breaks=seq(70, 100, by = 2), color="black")+
#   mytheme +
#   xlab("Coverage in 3' end (%)")+
#   ggtitle("Coverage of transcripts in 3' end")
# 
# # PLOT 22
# 
# p22 <- ggplot(data=perfect_match[!is.na(perfect_match$end5_perc_coverage_NA) & perfect_match$end5_perc_coverage_NA>=70,], aes(x=end5_perc_coverage_NA)) +
#   geom_histogram(fill=myPalette[3], breaks=seq(70, 100, by = 2), color="black") +
#   mytheme +
#   xlab("Coverage in 5' end (%)")+
#   ggtitle("Coverage of transcripts in 5' end")
#
# # PLOT 21.a
# 
# p21.a <- ggplot(data=perfect_match[perfect_match$diffToTTS>=0,], aes(x=diffToTTS)) +
#   geom_histogram(aes(y=..count../sum(..count..)), breaks = seq(0, 12000, by = 500), fill=myPalette[6], color="black")+
#   scale_y_continuous(labels = percent_format(), limits = c(0,1))+
#   mytheme +
#   xlab("# Nt lack to TTS")+
#   ggtitle("Lack of coverage in the 3'")  
# 
# p21.b <- ggplot(data=perfect_match[perfect_match$diffToTTS>=0 & perfect_match$diffToTTS<=1000,], aes(x=diffToTTS)) +
#   geom_histogram(aes(y=..count../sum(..count..)), breaks = seq(0, 500, by = 20), fill=myPalette[6], color="black")+
#   scale_y_continuous(labels = percent_format(), limits = c(0,1))+
#   mytheme +
#   xlab("# Nt lack to TTS")+
#   xlab("Distance to 3'end")+
#   ggtitle("Lack of coverage in the 3'")  
#
# p21.bb <- ggplot(data=perfect_match[which(perfect_match$diffToTTS<=250 & perfect_match$diffToTTS>=-250),], aes(x=diffToTTS)) +
#   geom_histogram(aes(y=..count../sum(..count..)), breaks = seq(-250, 250> dim(perfect_match)
#                                                                [1] 7775   19, by = 20), fill=myPalette[6], color="black")+
#   scale_y_continuous(labels = percent_format(), limits = c(0,1))+
#   mytheme +
#   xlab("# Nt lack to TTS")+
#   xlab("Distance to 3'end")+
#   ggtitle("Lack of coverage in the 3'")  
#
# p21.bbb <- ggplot(data=perfect_match[which(perfect_match$diffToTTS<=250 & perfect_match$diffToTTS>=-250),], aes(x=diffToTTS)) +
#   geom_histogram(breaks = seq(-250, 250, by = 20), fill=myPalette[6], color="black")+
#   mytheme +
#   xlab("# Nt lack to TTS")+
#   xlab("Distance to 3'end")+
#   ggtitle("Lack of coverage in the 3'")  
#
# ranges =  c(-10000,seq(-250, 250, by = 20),10000)
# p21.bbbb <- ggplot(data=perfect_match, aes(x=diffToTTS)) +
#   geom_histogram(breaks = ranges, fill=myPalette[6], color="black")+
#   mytheme +
#   xlab("# Nt lack to TTS")+
#   xlab("Distance to 3'end")+
#   ggtitle("Lack of coverage in the 3'")

roundUp = function(x) 10^ceiling(log10(x))
RoundUp = function(from,to) ceiling(from/to)*to

max = max(c(perfect_match$diffToTTS, perfect_match$diffToTSS)) 
perfect_match$rangeDiffTTS = cut(-(perfect_match$diffToTTS), breaks = c((max+1)*(-1),seq(-200, 200, by = 20),max+1))
perfect_match$rangeDiffTSS = cut(-(perfect_match$diffToTSS), breaks = c((max+1)*(-1),seq(-200, 200, by = 20),(max+1)))

max_height = max(c(max(table(perfect_match$rangeDiffTSS)),max(table(perfect_match$rangeDiffTTS))))
mm = roundUp(max_height)/10
max_height_rounded = RoundUp(max_height, mm)


p21.bbbb <- ggplot(data=perfect_match, aes(x=rangeDiffTTS)) +
  geom_bar(fill=myPalette[4], color="black", size=0.3)+
  scale_y_continuous(expand = c(0,0), limits = c(0,max_height_rounded))+
  mytheme +
  ylab("# FSM transcripts")+
  xlab("Distance to annotated TTS (bp)")+
  labs(     title="Distance distribution from sequenced to annotated TTS \n\n",
            subtitle="Negative values indicate that the sequenced TTS is downstream annotated TTS\n\n") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p21.bbbbb <- ggplot(data=perfect_match, aes(x=rangeDiffTTS)) +
  geom_bar(aes(y = (..count..)/sum(..count..)), fill=myPalette[4], color="black", size=0.3)+
  scale_y_continuous(labels = percent_format(), limits = c(0,1), expand = c(0,0))+
  mytheme +
  ylab("% FSM transcripts")+
  xlab("Distance to annotated TTS (bp)")+
  labs(     title="Distance distribution from sequenced to annotated TTS \n\n",
            subtitle="Negative values indicate that the sequenced TTS is downstream annotated TTS\n\n") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# p21.c <- ggplot(data=perfect_match[perfect_match$diffToTTS>=0 & perfect_match$diffToTTS<=1000,], aes(x=diffToTTS)) +
#   stat_ecdf(geom = "step", color=myPalette[6]) +
#   mytheme +
#   xlab("# Nt lack to TTS")+
#   ggtitle("Lack of coverage in the 3'")

# PLOT 22.a

# p22.a <- ggplot(data=perfect_match[perfect_match$diffToTSS>=0,], aes(x=diffToTSS)) +
#   geom_histogram(aes(y=..count../sum(..count..)), breaks = seq(0, 12000, by = 500), fill=myPalette[3], color="black")+
#   scale_y_continuous(labels = percent_format(), limits = c(0,1))+
#   mytheme +
#   xlab("# Nt lack to TSS")+
#   ggtitle("Lack of coverage in the 5'")
# 
# p22.b <- ggplot(data=perfect_match[perfect_match$diffToTSS>=0,], aes(x=diffToTSS)) +
#   geom_histogram(aes(y=..count../sum(..count..)), breaks = seq(0, 500, by = 20), fill=myPalette[3], color="black")+
#   scale_y_continuous(labels = percent_format(), limits = c(0,1))+
#   mytheme +
#   xlab("Distance to 5'end")+
#   ggtitle("Lack of coverage in the 5'")
# 
# p22.bb <- ggplot(data=perfect_match[which(perfect_match$diffToTSS<=250 & perfect_match$diffToTSS>=-250),], aes(x=diffToTSS)) +
#   geom_histogram(aes(y=..count../sum(..count..)), breaks = seq(-250, 250, by = 20), fill=myPalette[3], color="black")+
#   scale_y_continuous(labels = percent_format(), limits = c(0,1))+
#   mytheme +
#   xlab("Distance to 5'end")+
#   ggtitle("Lack of coverage in the 5'")
# 
# p22.bbb <- ggplot(data=perfect_match[which(perfect_match$diffToTSS<=250 & perfect_match$diffToTSS>=-250),], aes(x=diffToTSS)) +
#   geom_histogram(breaks = seq(-250, 250, by = 20), fill=myPalette[3], color="black")+
#   mytheme +
#   xlab("Distance to 5'end")+
#   ggtitle("Lack of coverage in the 5'")
#

p22.bbbbb <- ggplot(data=perfect_match, aes(x=rangeDiffTSS)) +
  geom_bar(aes(y = (..count..)/sum(..count..)), fill=myPalette[6], color="black", size=0.3)+
  scale_y_continuous(labels = percent_format(), limits = c(0,1), expand = c(0,0))+
  mytheme +
  ylab("% FSM transcripts")+
  xlab("Distance to annotated TSS (bp)")+
  labs(     title="Distance distribution from sequenced to annotated TSS \n\n",
            subtitle="Negative values indicate that the sequenced TSS is upsteam annotated TSS\n\n") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


p22.bbbb <- ggplot(data=perfect_match, aes(x=rangeDiffTSS)) +
  geom_bar(fill=myPalette[6], color="black", size=0.3)+
  scale_y_continuous(expand = c(0,0), limits = c(0,max_height_rounded))+
  mytheme +
  ylab("# FSM transcripts")+
  xlab("Distance to annotated TSS (bp)")+
  labs(     title="Distance distribution from sequenced to annotated TSS \n\n",
            subtitle="Negative values indicate that the sequenced TSS is upstream annotated TSS\n\n") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# p22.c <- ggplot(data=perfect_match[perfect_match$diffToTSS>=0 & perfect_match$diffToTSS<=1000,], aes(x=diffToTSS)) +
#   stat_ecdf(geom = "step", color=myPalette[3]) +
#   mytheme +
#   xlab("# Nt lack to TSS")+
#   ggtitle("Lack of coverage in the 5'")
# 
# 
# p22.q <- ggplot(data=perfect_match[perfect_match$diffToTSS>=0,], aes(x=diffToTSS)) +
#   geom_histogram(aes(y=..count../sum(..count..)), breaks = seq(0, 12000, by = 50), fill=myPalette[3], color="black")+
#   scale_y_continuous(labels = percent_format(), limits = c(0,1))+
#   mytheme +
#   xlab("# Nt lack to TSS")+
#   ggtitle("Lack of coverage in the 5'") +  geom_vline(xintercept = q_5)
# 
# 
# p21.q <- ggplot(data=perfect_match[perfect_match$diffToTTS>=0,], aes(x=diffToTTS)) +
#   geom_histogram(aes(y=..count../sum(..count..)), breaks = seq(0, 12000, by = 100), fill=myPalette[6], color="black")+
#   scale_y_continuous(labels = percent_format(), limits = c(0,1))+
#   mytheme +
#   xlab("# Nt lack to TTS")+
#   ggtitle("Lack of coverage in the 3'")   + geom_vline(xintercept = q_3)
# 
# # PLOT together

# 
# p21_22.a <- ggplot(data=to[to$diff>=0,], aes(x=diff, color=type)) +
#   stat_ecdf(geom = "step", size=0.8) +
#   mytheme +
#   ggtitle("Accumulative density") +
#   geom_vline(data = quantile95, aes(xintercept = diff, ylim = 1, color = type), linetype = 2, size=0.8, show.legend = FALSE) +
#   scale_color_manual(values = myPalette[c(3,6)]) +
#   theme(legend.box = 'horizontal',
#         legend.position = c(0.75, 0.75)) +
#   guides(fill = guide_legend(override.aes = list(linetype = 0))) +
#   annotate("text", mean(quantile95$diff), 1.05, label = "95th percentile")


# PLOT 23

p23 <- ggplot(data=junctionsData, aes(x=structuralCategory)) +
  geom_bar(position="fill", aes(y = (..count..)/sum(..count..), fill=canonical_known), color="black",  size=0.3, width = 0.7) +
  scale_y_continuous(labels = percent, expand = c(0,0)) +
  scale_fill_manual(values = myPalette[c(1,7,3,2)]) +
  ylab("% Splice junctions") +  
  mytheme +
  guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.3))+
  theme(legend.position="bottom", legend.title=element_blank())  +
  theme(axis.text.x = element_text(angle = 45)) + 
  theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12))+
  theme(axis.title.x=element_blank()) + 
  #theme(plot.margin = unit(c(1.5,1,0,1), "cm")) +
  ggtitle("Distribution of SJ type among structural classification\n\n\n") 

sqantiData$AllCanonical = factor(sqantiData$AllCanonical, 
                                       levels = c("canonical","non_canonical"), 
                                       ordered=TRUE)

p23.b <- ggplot(data=sqantiData, aes(x=structuralCategory)) +
  geom_bar(position="fill", aes(y = (..count..)/sum(..count..), fill=AllCanonical), color="black", size=0.3, width = 0.7) +
  scale_y_continuous(labels = percent, expand = c(0,0)) +
  scale_fill_manual(values = myPalette[c(1,7,3,2)]) +
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
p23.c <- ggplot(data=sqantiData, aes(x=structuralCategory, alpha=AllCanonical, fill=structuralCategory)) +
  geom_bar(position="fill", aes(y = (..count..)/sum(..count..)), color="black",  size=0.3,width = 0.7) +
  scale_fill_manual(values = myPalette, guide='none' )  +
  scale_y_continuous(labels = percent, expand = c(0,0)) +
  scale_alpha_manual(values=c(1,0.3)) +
  ylab("% Transcripts ") +  
  mytheme +
  theme(legend.position="bottom", axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45)) + 
  theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12))+
  theme(axis.title.x=element_blank()) + 
  #theme(plot.margin = unit(c(1.5,1,0,1), "cm")) +
  labs(   title="Distribution of transcripts by splice junction category\n\n\n",
          subtitle="Non canonical transcripts are those with at least one non-canonical junction\n\n") 
  

# PLOT 24

p24 <- ggplot(data=junctionsData, aes(x=transcriptCoord, fill = canonical_known)) +
  geom_density(alpha=0.7,  size=0.3) +
  scale_y_continuous(expand = c(0,0))  +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = myPalette[c(1,7,3,2)]) +
  geom_vline(xintercept = 200, linetype = "longdash", color="red") +
  mytheme +
  guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.7))+
  theme(legend.position="bottom" ) +
  labs(list(fill = "Junction Type", x= "Distance to TSS (bp)", 
            title = "Distribution of splice junctions distance to TSS\n\n") )


# PLOT 25

p25 <- ggplot(data=junctionsData, aes(x=TSSrange, fill=canonical_known)) +
  geom_bar(aes(y = (..count..)/sum(..count..)), color="black", size=0.3, width=0.7, position="fill") +
  scale_y_continuous(labels = percent, expand = c(0,0)) +
  scale_fill_manual(values = myPalette[c(1,7,3,2)]) +
  mytheme +
  theme(legend.position="bottom") +
  guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.7))+
  labs(list(fill = "Junction Type", x= "TSS distance range (bp)", y="% Junctions",
            title = "Splice junction distance to TSS across junction type\n\n\n") ) +
  theme(axis.title.x  = element_text(margin=margin(10,0,0,0), size=12))




# PLOT 26

p26 <- ggplot(data=junctionsData, aes(x=canonical_known, fill=TSSrange)) +
  geom_bar(aes(y = (..count..)/sum(..count..)), color="black", size=0.3, width=0.7, position="fill") +
  scale_y_continuous(labels = percent, expand = c(0,0)) +
  scale_fill_manual(values = myPalette) +
  ylab("% Junctions") +
  mytheme +
  theme(legend.position="bottom") +
  guides(fill = guide_legend(title = "TSS distance range (bp)")) +
  ggtitle( "Splice junction distance to TSS across junction type\n\n\n") +
  theme(axis.title.x=element_blank()) 
#theme(plot.margin = unit(c(1.5,1,0,1), "cm")) 

# PLOT 27

if (!all(is.na(sqantiData$MinCov))){
  
  junctionsData$isoExp = sqantiData[junctionsData$isoform, "isoExp"]
  junctionsData$junctionLabel = with(junctionsData, paste(chrom, strand,genomicStartCoord, genomicEndCoord, sep="_"))
  
  total = aggregate(cbind(totalCoverage,isoExp,transcriptCoord) ~ junctionLabel, data = junctionsData, 
                    FUN = function(x) c(mn = sum(x), n = min(x) ) )
  
  total$relCov = total$totalCoverage[,"n"] / total$isoExp[,"mn"]
  total$minTSS = total$transcriptCoord[,"n"]
  
  uniqJunc = unique(junctionsData[,c("junctionLabel", "canonical_known", "totalCoverage")])
  uniqJunc$notCov = uniqJunc$totalCoverage==0
  
  uniqueJunc_nonCov = as.data.frame(table(uniqJunc[uniqJunc$totalCoverage==0,"canonical_known"])/table(uniqJunc$canonical_known)*100)
  
  uniqJunc2 = merge(total, uniqJunc, by=1)
  uniqJunc2$TSSrange =cut(uniqJunc2$minTSS, breaks = c(0,40,80,120,160,200,10000000), labels = c("0-40", "41-80", "81-120", "121-160", "161-200",">200"))
  
  ### NO ENTIENDO!
  
  pn1 <- ggplot(data=uniqJunc2, aes(y=relCov,x=TSSrange,fill=canonical_known)) +
    geom_boxplot(outlier.size = 0.2) +
    scale_y_continuous(limits = c(0,1)) +
    scale_fill_manual(values = myPalette) +
    ylab("Relative coverage") + 
    xlab("# TSS distance range") +
    mytheme_bw +
    theme(legend.position="bottom", legend.title=element_blank())  +
    ggtitle( "Relative Coverage of junctions (unique junctions)\n\n\n")
    
  

  pn2 <-ggplot(data=uniqJunc2, aes(y=log(relCov+1),x=TSSrange,fill=canonical_known)) +
    geom_boxplot(outlier.size = 0.2) +
    scale_fill_manual(values = myPalette) +
    ylab("log(Relative coverage)") + 
    xlab("# TSS distance range") +
    mytheme_bw +
    theme(legend.position="bottom", legend.title=element_blank())  +
    ggtitle( "Relative Coverage of junctions (unique junctions) \n\n\n") 
  
  
  pn3 <-ggplot(data=uniqJunc2, aes(y=relCov,x=TSSrange,fill=canonical_known)) +
    geom_boxplot(outlier.size = 0.1) +
    scale_fill_manual(values = myPalette) +
    ylab("Relative coverage") + 
    xlab("# TSS distance range") +
    mytheme_bw +
    theme(legend.position="bottom", legend.title=element_blank())  +
    ggtitle( "Relative Coverage of junctions (unique junctions) \n\n\n") 
  
  pn4 <-ggplot(data=uniqJunc[uniqJunc$notCov==TRUE,], aes(x=canonical_known,fill=canonical_known)) +
    geom_bar(color="black", size=0.3, width=0.7) +
    scale_fill_manual(values = myPalette[c(1,7,3,2)]) +
    ylab("# Splice Junctions without coverage") + 
    xlab("Junction Type") +
    mytheme +
    guides(fill=FALSE) +
    ggtitle( "Splice junctions without short-read coverage (unique junctions) \n\n\n") 

  
  pn5 <- ggplot(data=uniqueJunc_nonCov, aes(x=Var1, y=Freq, fill=Var1)) +
    geom_bar(position=position_dodge(), stat="identity", color="black", size=0.3, width=0.7) +
    guides(fill=FALSE) +
    scale_y_continuous( expand = c(0,0)) +
    scale_fill_manual(values = myPalette[c(1,7,3,2)]) +
    ylab("% Splice Junctions without coverage") + 
    xlab("Junction Type") +
    mytheme +
    guides(fill=FALSE) +
    ggtitle( "Splice junctions without short-read coverage (unique junctions) \n\n\n") 
  
  # same but using all the junctions
  
  junctionsData$isoExp = sqantiData[junctionsData$isoform, "isoExp"]
  junctionsData$junctionLabel = with(junctionsData, paste(chrom, strand,genomicStartCoord, genomicEndCoord, sep="_"))
  
  # calculate sum of mean expression of transcripts where each junction
  sumExpPerJunc = tapply(junctionsData$isoExp, junctionsData$junctionLabel, sum)  #71387 uniq junction
  
  junctionsData$sumIsoExp = sumExpPerJunc[junctionsData$junctionLabel]
  
  junctionsData$relCov = junctionsData$totalCoverage / junctionsData$sumIsoExp
  
  junctionsData$TSSrange =cut(junctionsData$transcriptCoord, breaks = c(0,20,40,60,80,100,120,140,160,180,200,10000000), labels = c("0-20", "21-40","41-80","61-80", "81-100","101-120", "121-140","141-160", "161-180", "181-200", ">200"))
  
  
  pn1.2 <-ggplot(data=junctionsData[junctionsData$relCov<1,], aes(y=relCov,x=TSSrange,fill=canonical_known)) +
    geom_boxplot(outlier.size = 0.2, size=0.3) +
    scale_fill_manual(values = myPalette[c(1,7,3,2)]) +
    ylab("Relative coverage") + 
    xlab("# TSS distance range") +
    mytheme_bw +
    theme(legend.position="bottom", legend.title=element_blank())  +
    ggtitle( "Relative Coverage of junctions\n\n\n") +
    theme(axis.text.x = element_text(angle = 45,margin=margin(15,0,0,0), size=12)) 
    
  
  
  pn2.2 <-ggplot(data=junctionsData, aes(y=log(relCov+1),x=TSSrange,fill=canonical_known)) +
    geom_boxplot(outlier.size = 0.1) +
    scale_fill_manual(values = myPalette) +
    ylab("log(Relative coverage)") + 
    xlab("# TSS distance range") +
    mytheme_bw +
    theme(legend.position="bottom", legend.title=element_blank())  +
    ggtitle( "Relative Coverage of junctions\n\n\n") 
  
  
  pn3.2 <-ggplot(data=junctionsData, aes(y=relCov,x=TSSrange,fill=canonical_known)) +
    geom_boxplot(outlier.size = 0.1) +
    scale_fill_manual(values = myPalette) +
    ylab("Relative coverage") + 
    xlab("# TSS distance range") +
    mytheme_bw +
    theme(legend.position="bottom", legend.title=element_blank())  +
    ggtitle( "Relative Coverage of junctions\n\n\n") 
  
}

# PLOT 28

p28 <- ggplot(data=attribute.df2[attribute.df2$Var1%in%c("FSM", "NNC", "NIC"),], aes(x=Var1, y=perc, fill=Var1, alpha= Var2)) +
  geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = myPalette[c(1,3,4)], guide='none' ) +
  scale_alpha_manual(values=c(1,0.4,0.1), name = "QC Attributes") +
  ylab("% Transcripts") + 
  xlab("") +
  mytheme +
  theme(legend.position="bottom", axis.title.x = element_blank()) +
  ggtitle( "Quality control attributes across structural categories\n\n" ) 


p28.a <- ggplot(data=attribute.df2[attribute.df2$Var1%in%c("FSM", "NNC", "NIC"),], aes(x=Var1, y=perc, fill= Var2)) +
  geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = myPalette[9:11]) +
  ylab("% Transcripts") + 
  xlab("") +
  mytheme +
  theme(legend.position="bottom", axis.title.x = element_blank()) +
  ggtitle( "Quality control attributes across structural categories\n\n" ) +
  guides(fill = guide_legend(title = "QC Attributes") )



attribute.df.2.2$structuralCategory = unlist(lapply(strsplit(as.character(attribute.df.2.2$Var1), split=":"), function(x) x[[2]]))


p28.c <- ggplot(data=attribute.df.2.2[attribute.df.2.2$structuralCategory%in%c("FSM","NNC"),], aes(x=Var1, y=perc, fill=Var2)) +
  geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = myPalette[9:11]) +
  ylab("% Transcripts") + 
  xlab("") +
  mytheme +
  theme(legend.position="bottom", axis.title.x = element_blank()) +
labs ( title= "Quality control attributes across structural categories\n\n",
         subtitle="Categories are divided into NNC/FSM containing genes") + 
  guides(fill = guide_legend(title = "QC Attributes") ) +
  scale_x_discrete(breaks=c("Genes expressing\nboth NNC and\n FSM isoforms:FSM", 
                          "Genes expressing\nboth NNC and\n FSM isoforms:NNC",
                          "Genes expressing\n only FSM isoforms:FSM",
                          "Genes expressing\n only NNC isoforms:NNC"),
                 labels=c("FSM\nNNC/FSM genes", 
                          "NNC\nNNC/FSM genes",
                          "FSM\nFSM genes",
                          "NNC\nNNC genes")) +
  theme(axis.text.x = element_text(size=10))
        



p28.d <- ggplot(data=attribute.df.2.2, aes(x=Var1, y=perc, fill=Var2)) +
  geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values =  myPalette[9:11]) +
  ylab("% Transcripts") + 
  xlab("") +
  mytheme +
  theme(legend.position="bottom", axis.title.x = element_blank()) +
  ggtitle( "Quality control attributes across structural categories\n\n" ) +
  guides(fill = guide_legend(title = "QC Attributes") ) +
  scale_x_discrete(breaks=c("Genes expressing\nboth NNC and\n FSM isoforms:FSM", 
                            "Genes expressing\nboth NNC and\n FSM isoforms:Genic\nGenomic",
                            "Genes expressing\nboth NNC and\n FSM isoforms:ISM",
                            "Genes expressing\nboth NNC and\n FSM isoforms:NIC",
                            "Genes expressing\nboth NNC and\n FSM isoforms:NNC",
                            "Genes expressing\n only FSM isoforms:FSM",
                            "Genes expressing\n only NNC isoforms:NNC"),
                   labels=c("FSM\nNNC/FSM genes", 
                            "Genic Genomic\nNNC/FSM genes",
                            "ISM\nNNC/FSM genes",
                            "NIC\nNNC/FSM genes",
                            "NNC\nNNC/FSM genes",
                            "FSM\nFSM genes",
                            "NNC\nNNC genes"))










# PLOT 29

b = table(junctionsData$canonical, junctionsData$junctionCategory)
a = table(junctionsData[junctionsData$RTS_junction=="TRUE","canonical"], junctionsData[junctionsData$RTS_junction=="TRUE","junctionCategory"])
rts_junction = as.data.frame(a/b*100)
rts_junction$`Junction Type` = with(rts_junction, paste(Var2,Var1, sep="\n"))
rts_junction$`Junction Type` = factor(rts_junction$`Junction Type`, levels=c("known\ncanonical", "known\nnon_canonical","novel\ncanonical","novel\nnon_canonical"),
                                       labels=c("Known\ncanonical ", "Known\nNon-canonical ", "Novel\ncanonical ", "Novel\nNon-canonical "), order=T) 

max_height_rts = max(rts_junction$Freq)
mm_rts = roundUp(max_height_rts)/10
max_height_rounded_rts = RoundUp(max_height_rts, mm_rts)

p29 <- ggplot(data=rts_junction, aes(x=`Junction Type`, y=Freq, fill=`Junction Type`)) +
  geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
  scale_fill_manual(values = myPalette[c(1,7,3,2)]) +
  ylab("% RT-switching junctions") +
  ggtitle( "RT-switching by splice junction category \n\n" ) +
  mytheme +
  guides(fill=FALSE) +
  scale_y_continuous(expand = c(0,0), limits = c(0,max_height_rounded_rts)) +
  theme(axis.text.x = element_text(size=11))

  



uniqJuncRTS =unique(junctionsData[,c("junctionLabel","canonical_known", "RTS_junction")])
c = table(uniqJuncRTS$canonical_known)
d = table(uniqJuncRTS[uniqJuncRTS$RTS_junction==TRUE, "canonical_known"])
rts_junction_uniq = as.data.frame(d/c*100)

max_height_rts_u = max(rts_junction_uniq$Freq)
mm_rts_u = roundUp(max_height_rts_u)/10
max_height_rounded_rts_u = RoundUp(max_height_rts_u, mm_rts_u)

p29.a <- ggplot(data=rts_junction_uniq, aes(x=Var1, y=Freq, fill=Var1)) +
  geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
  scale_fill_manual(values = myPalette[c(1,7,3,2)]) +
  ylab("% RT-switching junctions") +
  ggtitle( "RT-switching by splice junction category (unique junctions)\n\n" ) +
  mytheme +
  guides(fill=FALSE) +
  scale_y_continuous(expand = c(0,0), limits = c(0,max_height_rounded_rts_u)) +
  theme(axis.text.x = element_text(size=11))




###** Output plots

pdf(file=report.file, width = 6.5, height = 6.5)


#cover
grid.newpage()
cover <- textGrob("SQANTI report",
                  gp=gpar( fontface="italic", fontsize=40, col="orangered"))
grid.draw(cover)

# t1
freqCat = as.data.frame(table(sqantiData$structuralCategory))
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
freqCat = as.data.frame(table(uniqJunc$canonical_known))
freqCat$Var1 = gsub(" ", "", freqCat$Var1)
freqCat$Var1 = gsub("\n", " ", freqCat$Var1)
table2 <- tableGrob(freqCat, rows = NULL, cols = c("category","# SJ"))
title2 <- textGrob("SJ classification", gp=gpar(fontface="italic", fontsize=17), vjust = -5)
gt3 <- gTree(children=gList(table2, title2))

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
#p0.2
p0
p6
p7
p10
if (!all(is.na(sqantiData$FL))){
  p11
}

#2. general parameters by structual categories

s <- textGrob("Structrual Isoform characterization\nbased on splice junctions", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
grid.arrange(s)
p1
p4
p5
p8
if (!all(is.na(sqantiData$FL))){
  print(p9)
}
p2
p3
p13
p13.c
p12



#3. splice junction

s <- textGrob("Splice junction characterization", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
grid.arrange(s)
p23.c
p23.b
p23
p24
p25
p26
if (!all(is.na(sqantiData$MinCov))){
  print(pn1.2)
  print(pn4)
  print(pn5)
}
p29
p29.a

#4. full-lengthness

s <- textGrob("Full-lengthness characterization of isoforms", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
grid.arrange(s)
p21.bbbb
p21.bbbbb
p22.bbbb
p22.bbbbb



s <- textGrob("Quality control attributes", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
grid.arrange(s)
#p28a
p28.a
p28.c


dev.off()



