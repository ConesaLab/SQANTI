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

#class.file = "/home/ldelafuente/squanti/sqanti_test/all.good.5merge.collapsed.longest_rep_classification.txt"
#junc.file = "/home/ldelafuente/squanti/sqanti_test//all.good.5merge.collapsed.longest_rep_junctions.txt"

#********************** Packages (installed if not found)

packages = c("ggplot2", "scales", "reshape")

for (package in packages){
  if (!package %in% installed.packages()) install.packages(package)
}


library(ggplot2)
library(scales)
library(reshape)


#********************** Reading Information

########## Classification information

sqantiData = read.table(file=class.file, header=T, as.is=T, sep="\t")
rownames(sqantiData) = sqantiData$isoform
xaxislabelsF1 = c("FSM", "ISM", "NIC", "NNC", "Intergenic", "Genic\nIntron", "Genic\nGenomic",  "Antisense", "Fusion")
sqantiData$structuralCategory = factor(sqantiData$structuralCategory, 
                                       labels = xaxislabelsF1, 
                                       levels = c("full-splice_match","incomplete-splice_match","novel_in_catalog","novel_not_in_catalog","intergenic","genic_intron","genic","antisense","fusion"), 
                                       ordered=TRUE)

########### Junction information

junctionsData = read.table(file=junc.file, header=T, as.is=T, sep="\t")

########### Managing information 

sqantiData[grep("novelGene",  sqantiData$associatedGene), "novelGene"] <- "Novel Genes"
sqantiData[is.na(sqantiData$novelGene), "novelGene"] <- "Annotated Genes"
#sqantiData[-grep("novelGene",  sqantiData$associatedGene), "novelGene"] <-"Annotated\nGenes"
sqantiData[which(sqantiData$exons>1), "exonCat"] <- "MultiExon"
sqantiData[-which(sqantiData$exons>1), "exonCat"] <- "MonoExon"

sqantiData[which(sqantiData$MinCov==0), "Not coverage SJ"] <- "Not coverage SJ"
a = sqantiData[,c("AllCanonical", "Not coverage SJ", "RTS_stage", "isoform", "structuralCategory")]
#a = sqantiData[,c("AllCanonical", "RTS_stage", "isoform", "structuralCategory")]
a$RTS_stage = as.factor(a$RTS_stage)
a$`Not coverage SJ` = as.factor(a$`Not coverage SJ`)
b = melt(a, id.vars = c("isoform", "structuralCategory"), na.rm = T)
b = b[-which(b$value=="FALSE" | b$value=="canonical"),]
b$value <- factor(b$value)
b$value = factor(b$value, 
                 labels = c("Non-canonical SJ","RT-swiching", "Not coverage SJ"), 
                 levels = c("non_canonical", "TRUE", "Not coverage SJ"), 
                 ordered=TRUE)
attribute.df = table(b$structuralCategory,b$value)
attribute.df2 = merge(as.data.frame(attribute.df), as.data.frame(table(sqantiData$structuralCategory)), by = 1)
attribute.df2$perc = attribute.df2$Freq.x / attribute.df2$Freq.y * 100

junctionsData$canonical_known = with(junctionsData, paste(junctionCategory,canonical,"SJ", sep="_"))
junctionsData$canonical_known=as.factor(junctionsData$canonical_known)
junctionsData$canonical_known = factor(junctionsData$canonical_known, levels=c("known_canonical_SJ", "known_non_canonical_SJ", "novel_canonical_SJ", "novel_non_canonical_SJ"),
                                       labels=c("Known\ncanonical ", "Known\nNon-canonical ", "Novel\ncanonical ", "Novel\nnon-canonical "), order=T) 
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

FSM_genes = unique(sqantiData[sqantiData$FSM_class=="A" & sqantiData$structuralCategory=="FSM","associatedGene"])  # unique isoform is a FSM
NNC_genes = unique(sqantiData[sqantiData$FSM_class=="A" & sqantiData$structuralCategory=="NNC","associatedGene"])  # unique isoform is a FSM #1520 instead of 585
FSMandNNCgenes = unique(sqantiData[sqantiData$FSM_class=="C" & sqantiData$structuralCategory=="NNC","associatedGene"]) # a FSM and any NNC isoform

gene[gene$associatedGene %in% FSMandNNCgenes, "class"] <- "Genes expressing\nboth NNC and\n FSM isoforms"
gene[gene$associatedGene %in% NNC_genes, "class"] <- "Genes expressing\n only NNC isoforms"
gene[gene$associatedGene %in% FSM_genes, "class"] <- "Genes expressing\n only FSM isoforms"

sqantiData[sqantiData$associatedGene %in% FSMandNNCgenes, "class"] <- "Genes expressing\nboth NNC and\n FSM isoforms"
sqantiData[sqantiData$associatedGene %in% NNC_genes, "class"] <- "Genes expressing\n only NNC isoforms"
sqantiData[sqantiData$associatedGene %in% FSM_genes, "class"] <- "Genes expressing\n only FSM isoforms"


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


# NEW: including also longer ends --> diff = 0
#perfect_match[which(perfect_match$diffToTSS<0),"diffToTSS"] <- 0
#perfect_match[which(perfect_match$diffToTTS<0),"diffToTTS"] <- 0


#perfect_match$range3 =cut(perfect_match$end3_perc_coverage_NA, breaks = seq(0,100,by = 5))

q_5 = quantile(perfect_match[perfect_match$diffToTSS>=0,"diffToTSS"], probs = seq(0,1,0.1))
q_3 = quantile(perfect_match[perfect_match$diffToTTS>=0,"diffToTTS"], probs = seq(0,1,0.1))
to = data.frame(diff = c(perfect_match$diffToTSS, perfect_match$diffToTTS), type = rep(c("Nt to TSS", "Nt to TTS"), each=nrow(perfect_match)))
to = to[to$diff>=0,]
quantile95 = aggregate(to[1], to[c(2)], function(x) quantile(x,0.95))


#*** Generating plots
# library(RColorBrewer)
# colo = NULL
# for (i in rn){
#   colo =c(colo,brewer.pal(9, name=i)[5])
# }
# 
# ej = colo[c(-5,-9,-7,-11,-3,-18)]
# ej2 = c(ej, "#E69F00","coral2")
# ej3 = ej2[c(-7,-8,-9)]

# myPalette1 <-  c("#66C2A5", "#FC8D62", "#6BAED6", "#FEE08B" ,"#7FC97F", "#FBB4AE","#3288BD", "#F46D43",  "#B3CDE3", "#CCEBC5" ,"#DECBE4")
# myPalette2 = c("#66C2A4","#FC8D59","#6BAED6","#67A9CF","#41B6C4","coral2","#78C679","#969696","#FE9929" ,"#7BCCC4" ,"goldenrod1" )
# myPalette3 = c("#78C679","#FC8D59","#6BAED6","#66C2A4","#41B6C4","coral2","#67A9CF","#969696","#FE9929" ,"#7BCCC4" ,"goldenrod1" )
# myPalette4 = c("#66C2A4","#FC8D59","#6BAED6","#78C679","coral2","#41B6C4","#969696","#FE9929", "#7BCCC4" ,"goldenrod1", "#67A9CF" )
# myPalette5 = c("#66C2A4","#FC8D59","#6BAED6","#78C679","coral2","#41B6C4","#969696", "goldenrod1", "#7BCCC4" , "#FE9929","#67A9CF" )
# myPalette6 = c("#66C2A4","#FC8D59","#6BAED6","coral2","#78C679","#969696", "goldenrod1", "#7BCCC4" , "#FE9929", "#41B6C4", "#67A9CF")
# myPalette = c("#66C2A4","#FC8D59","#6BAED6","coral2","#78C679","#7BCCC4", "goldenrod1",  "#969696", "#FE9929", "#41B6C4", "#67A9CF")
# #myPalette = c("#66C2A4","#FC8D59","#6BAED6","coral2","#78C679","#7BCCC4", "goldenrod1",  "#969696",  "#41B6C4", "#67A9CF")
# myPalette = c("#66C2A4","#FC8D59","#6BAED6","goldenrod1","#78C679","#7BCCC4", "coral2",  "#969696", "#FE9929", "#41B6C4", "#67A9CF")

myPalette = c("#66C2A4","#FC8D59","#6BAED6","goldenrod1","#78C679","coral2","#969696", "#7BCCC4" , "#FE9929", "#41B6C4", "#67A9CF")

mytheme <- theme_classic(base_family = "Helvetica") +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4)) +
  theme(axis.title.x = element_text(size=14, margin=margin(5,0,0,0)),
        axis.text.x  = element_text(margin=margin(7,0,0,0), size=13),
        axis.title.y = element_text(size=14,  margin=margin(0,15,0,0)),
        axis.text.y  = element_text(vjust=0.5, size=13) ) +
  theme(legend.text = element_text(size = 10), legend.title = element_text(size=11), legend.key.size = unit(0.5, "cm")) +
  theme(plot.title = element_text(lineheight=.4, size=15, margin=margin(-20,0,0,0))) +
  theme(legend.justification=c(1,1), legend.position=c(1,1)) +
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
  theme(plot.title = element_text(lineheight=.4, size=15, margin=margin(-20,0,0,0))) +
  theme(legend.justification=c(1,1), legend.position=c(1,1)) +
  theme(plot.margin = unit(c(2.5,1,1,1), "cm")) 

#myPalettes = list(myPalette1, myPalette2, myPalette3, myPalette4, myPalette5, myPalette6)

#for (i in 1:length(myPalettes)){
  # print(i)
  # myPalette = myPalettes[[i]]

# print(ggplot(data=sqantiData, aes(x=structuralCategory, y=length, fill=structuralCategory)) +
#     geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
#     xlab("Structural Classification") +  
#     ylab("Transcript Length (bp)") +
#     scale_fill_manual(values = myPalette) +
#     guides(fill=FALSE) +
#     mytheme  +
#     ggtitle("Transcript length distribution\n\n" ))
#   
# }

# PLOT 1

legendLabelF1 = tools::toTitleCase(gsub("_", " ", levels(as.factor(sqantiData$coding))))

p1 <- ggplot(data=sqantiData, aes(x=structuralCategory)) +
  geom_bar(aes(y = (..count..)/sum(..count..), fill=coding), color="black", size=0.3, width=0.7) +
  scale_y_continuous(labels = percent, expand = c(0,0)) +
  geom_text(aes(y = ((..count..)/sum(..count..)), label = scales::percent((..count..)/sum(..count..))), stat = "count", vjust = -0.25)  +
  scale_fill_manual(name = "Coding prediction", 
                    labels = legendLabelF1, 
                    values = myPalette) +
  xlab("") + 
  ylab("% Transcripts") +
  mytheme + 
  geom_blank(aes(y=1.05*((..count..)/sum(..count..))), stat = "count") +
  theme(axis.text.x = element_text(angle = 45)) +
  ggtitle("Isoform distribution across Structural Categories\n\n" ) +
  theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12)) 
  #theme(plot.margin = unit(c(1.5,1,0,1), "cm")) 




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
  ggtitle("Transcript length distribution by Structural Classification\n\n" ) +  
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
  ggtitle("Exon number distribution by Structural Classification\n\n" ) +
  theme(axis.title.x=element_blank()) 
  #theme(plot.margin = unit(c(1.5,1,0,1), "cm")) 


# PLOT 6

p6 <- ggplot(data=sqantiData, aes(x=novelGene)) +
  geom_bar(position="fill",aes(y = (..count..)/sum(..count..), fill=exonCat), color="black", size=0.3, width=0.5) +
  #geom_text(aes(y = (..count..), label=(..count..)), stat = "count", vjust = -0.25)  +
  scale_y_continuous(labels = percent, expand = c(0,0)) +
  scale_fill_manual(name = "Transcript type", 
                    values = myPalette) +
  ylab("% Transcripts ") +  
  mytheme +
  theme(axis.title.x=element_blank()) + 
  theme(legend.position="bottom") +   
  ggtitle("Distribution of mono and multi exon transcripts\n\n by type of Gene Annotation\n\n" ) 



# PLOT 7

p7 <- ggplot(data=gene, aes(x=novelGene)) +
  geom_bar(position="fill", aes(y = (..count..)/sum(..count..), fill=range), color="black", size=0.3, width=0.5) +
  scale_y_continuous(labels = percent, expand = c(0,0)) +
  scale_fill_manual(name = "# Isoforms Per Gene",
                    values = myPalette) +
  ylab("% Genes ") +  
  xlab("Gene Type") +
  mytheme +
  theme(axis.title.x=element_blank()) + 
  theme(legend.position="bottom") +   
  guides(fill = guide_legend(keywidth = 0.9, keyheight = 0.9)) +
  ggtitle("Distribution of number of isoforms by type of Gene Annotation\n\n\n\n" ) 
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
  ggtitle("Transcript Expression by Structural Category\n\n" )


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
    ggtitle("Number of FL reads per transcript by Structural Category\n\n" )

}

# PLOT 10 

p10 <- ggplot(data=gene, aes(x=novelGene, y=log2(geneExp+1), fill=novelGene)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  xlab("Structural Classification") +  
  ylab("log2( # Short reads + 1)") +
  scale_fill_manual(values = myPalette) +
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
    scale_fill_manual(values = myPalette) +
    guides(fill=FALSE) +
    mytheme_bw +
    theme(axis.title.x=element_blank()) + 
    #theme(plot.margin = unit(c(1.5,1,0.5,1), "cm")) +
    ggtitle("Number of FL reads per Gene by type of Gene Annotation\n\n" )

}

# PLOT 12

p12 <- ggplot(data=gene, aes(x=FSM_class2, y=log2(geneExp+1), fill=FSM_class2)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  xlab("") +  
  ylab("log2(# Short reads per gene + 1)") +
  scale_fill_manual(values = myPalette) +
  guides(fill=FALSE) +
  mytheme_bw +
  theme(axis.title.x=element_blank()) + 
  
  ggtitle("Transcript expression levels across A/B/C groups\n\n" )


# PLOT 13

p13 <- ggplot(data=gene, aes(x=class, y=log2(geneExp+1), fill=class)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  ylab("log2( # Short reads per gene + 1)") +
  theme(axis.title.x=element_blank()) + 
  #theme(plot.margin = unit(c(1.5,1,0.5,1), "cm")) +
  scale_fill_manual(values = myPalette) +
  guides(fill=FALSE) + 
  mytheme_bw +
  theme(axis.title.x=element_blank()) + 
  ggtitle("Gene Expression level in NNC and FSM containing genes\n\n" )


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

p15 <- ggplot(data=sqantiData, aes(x=class, y=log2(isoExp+1), fill=structuralCategory)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  ylab("log2( # Short reads per transcript + 1)") +
  #theme(plot.margin = unit(c(1.5,1,0.5,1), "cm")) +
  scale_fill_manual(values = myPalette) +
  mytheme +
  theme(axis.title.x=element_blank()) + 
  ggtitle("Transcript expression levels for NNC and FSM containing genes\n\n" ) +
  guides(fill = guide_legend(title = "Structural Classification", keywidth = 0.9, keyheight = 0.9)) +
  theme(legend.position="bottom")

  
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
# # PLOT 20
# 
# p20 <- ggplot(data=perfect_match[!is.na(perfect_match$end5_perc_coverage_NA) & perfect_match$end5_perc_coverage_NA>=70,], aes(x=end5_perc_coverage_NA)) +
#   geom_histogram(aes(y=..density..),fill="white", breaks=seq(70, 100, by = 2), color="black") +
#   geom_density(alpha=.6, fill=myPalette[3]) +
#   mytheme +
#   xlab("Coverage in 5' end (%)")+
#   ylab("Density") +
#   ggtitle("Coverage of transcripts in 5' end")

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

  
# p21.bb <- ggplot(data=perfect_match[which(perfect_match$diffToTTS<=250 & perfect_match$diffToTTS>=-250),], aes(x=diffToTTS)) +
#   geom_histogram(aes(y=..count../sum(..count..)), breaks = seq(-250, 250> dim(perfect_match)
#                                                                [1] 7775   19, by = 20), fill=myPalette[6], color="black")+
#   scale_y_continuous(labels = percent_format(), limits = c(0,1))+
#   mytheme +
#   xlab("# Nt lack to TTS")+
#   xlab("Distance to 3'end")+
#   ggtitle("Lack of coverage in the 3'")  

  
# p21.bbb <- ggplot(data=perfect_match[which(perfect_match$diffToTTS<=250 & perfect_match$diffToTTS>=-250),], aes(x=diffToTTS)) +
#   geom_histogram(breaks = seq(-250, 250, by = 20), fill=myPalette[6], color="black")+
#   mytheme +
#   xlab("# Nt lack to TTS")+
#   xlab("Distance to 3'end")+
#   ggtitle("Lack of coverage in the 3'")  

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
  geom_bar(fill=myPalette[6], color="black", size=0.3)+
  scale_y_continuous(expand = c(0,0), limits = c(0,max_height_rounded))+
  mytheme +
  ylab("# FSM transcripts")+
  xlab("Distance to annotated TTS (bp)")+
  ggtitle("Distance distribution from sequenced to annotated TTS\n\n")  +
  theme(axis.text.x =element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank())

p21.bbbbb <- ggplot(data=perfect_match, aes(x=rangeDiffTTS)) +
  geom_bar(aes(y = (..count..)/sum(..count..)), fill=myPalette[6], color="black", size=0.3)+
  scale_y_continuous(labels = percent_format(), limits = c(0,1), expand = c(0,0))+
  mytheme +
  ylab("% FSM transcripts")+
  xlab("Distance to annotated TTS (bp)")+
  ggtitle("Distance distribution from sequenced to annotated TTS\n\n") +
  theme(axis.text.x =element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank())



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



p22.bbbb <- ggplot(data=perfect_match, aes(x=rangeDiffTSS)) +
  geom_bar(fill=myPalette[3], color="black", size=0.3)+
  scale_y_continuous(expand = c(0,0), limits = c(0,max_height_rounded))+
  mytheme +
  ylab("# FSM transcripts")+
  xlab("Distance to annotated TSS (bp)")+
  ggtitle("Distance distribution from sequenced to annotated TSS\n\n")  +
  theme(axis.text.x =element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank())

p22.bbbbb <- ggplot(data=perfect_match, aes(x=rangeDiffTSS)) +
  geom_bar(aes(y = (..count..)/sum(..count..)), fill=myPalette[3], color="black", size=0.3)+
  scale_y_continuous(labels = percent_format(), limits = c(0,1), expand = c(0,0))+
  mytheme +
  ylab("% FSM transcripts")+
  xlab("Distance to annotated TSS (bp)")+
  ggtitle("Distance distribution from sequenced to annotated TSS\n\n") +
  theme(axis.text.x =element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank())

#   
# scale_x_discrete(breaks=c(unique(perfect_match$rangeDiffTSS)),
#                    labels=c(rep("",22)))
#                    
#                       labels=c("-10000","-200",rep("",18),"+200","+10000")) 


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
  geom_bar(position="fill", aes(y = (..count..)/sum(..count..), fill=canonical_known), color="black", width = 0.7) +
  scale_y_continuous(labels = percent, expand = c(0,0)) +
  scale_fill_manual(values = myPalette) +
  ylab("% Transcripts ") +  
  mytheme +
  guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.3))+
  theme(legend.position="bottom", legend.title=element_blank())  +
  theme(axis.text.x = element_text(angle = 45)) + 
  theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12))+
  ggtitle("Exon number distribution by Structural Classification\n\n" ) +
  theme(axis.title.x=element_blank()) + 
  #theme(plot.margin = unit(c(1.5,1,0,1), "cm")) +
  ggtitle("Distribution of SJ type among Structural Classification\n\n\n") 
  


# PLOT 24

p24 <- ggplot(data=junctionsData, aes(x=transcriptCoord, fill = canonical_known)) +
  geom_density(alpha=0.7) +
  scale_y_continuous(expand = c(0,0))  +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = myPalette) +
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
  scale_fill_manual(values = myPalette) +
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
  guides(fill = guide_legend(title = "TSS distance range (bp)", title.position = "left")) +
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
  
  uniqJunc = unique(junctionsData[,c("junctionLabel", "canonical_known")])
  
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
    geom_boxplot(outlier.size = 0.2) +
    scale_fill_manual(values = myPalette) +
    ylab("Relative coverage") + 
    xlab("# TSS distance range") +
    mytheme_bw +
    theme(legend.position="bottom", legend.title=element_blank())  +
    ggtitle( "Relative Coverage of junctions (unique junctions) \n\n\n") 
    
  
  
  # same but using all the junctions
  
  
  
  junctionsData$isoExp = sqantiData[junctionsData$isoform, "isoExp"]
  junctionsData$junctionLabel = with(junctionsData, paste(chrom, strand,genomicStartCoord, genomicEndCoord, sep="_"))
  
  junctionsData$relCov = junctionsData$totalCoverage / junctionsData$isoExp
  
  junctionsData$TSSrange =cut(junctionsData$transcriptCoord, breaks = c(0,40,80,120,160,200,10000000), labels = c("0-40", "41-80", "81-120", "121-160", "161-200",">200"))
  
  
  pn1.2 <-ggplot(data=junctionsData, aes(y=relCov,x=TSSrange,fill=canonical_known)) +
    geom_boxplot(outlier.size = 0.1) +
    scale_y_continuous(limits = c(0,1)) +
    scale_fill_manual(values = myPalette) +
    ylab("Relative coverage") + 
    xlab("# TSS distance range") +
    mytheme_bw +
    theme(legend.position="bottom", legend.title=element_blank())  +
    ggtitle( "Relative Coverage of junctions\n\n\n") 
    
  
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

p28 <- ggplot(data=attribute.df2[attribute.df2$Var1%in%c("FSM", "NNC", "NIC"),], aes(x=Var1, y=perc, fill= Var2)) +
  geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = myPalette) +
  ylab("% Transcripts in each class ") + 
  xlab("") +
  mytheme +
  theme(legend.position="bottom", axis.title.x = element_blank()) +
  ggtitle( "Quality Control Attributes accros Structural Categories\n\n" ) +
  guides(fill = guide_legend(title = "QC Attributes") )



### Output plots

pdf(file=report.file, width = 6.5, height = 6.5)

#pdf(file=paste("Report_17OCT",i,sep="_"),  width = 7, height = 7)
# 
# print(p1)
# print(p2)
# print(p3)
# print(p4)
# print(p5)
# print(p6)
# print(p7)
# print(p8)
# print(p9)
# print(p10)
# print(p11)
# print(p12)
# print(p13)
# print(p15)
# print(p16)
# print(p17)
# print(p18)
# print(p19)
# print(p20)
# print(p21)
# print(p22)
# print(p23)
# print(p24)
# print(p25)
# print(p26)
# print(p28)
# 
# dev.off()
# 
# print("SQANTI report finished successfully.")
# 
# 
# 
# 
# }


#pdf(file="Report_20OCT.pdf",  width = 7, height = 7)


p1
p28
p2
p3
p4
p5
p6
p7
p8
if (!all(is.na(sqantiData$FL))){
  print(p9)
  print(p11)
}
p10
p12
p13
p15
p21.bbbb
p22.bbbb
p21.bbbbb
p22.bbbbb
p23
p24
p25
p26
if (!all(is.na(sqantiData$MinCov))){
  print(pn1)
  print(pn2)
  print(pn3)
  print(pn1.2)
  print(pn2.2)
  print(pn3.2)
}

dev.off()


print("SQANTI QC succesfully done!")
