# Cecile Pereira & Lorena de la Fuente
# cecilepereira@ufl.edu

### package requirements

list.of.packages <- c("rpart", "ROCR", "caret", "optparse", "ggplot2", "lattice", "foreach", "e1071","randomForest", "partykit", "ipred", "rpart.plot", "doMC", "nnet", "ROSE", "pROC", "MLmetrics")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

rm(list=ls())
library(rpart)           #Rpart package
library(partykit)        #"Elegant" tree design package
library(ipred)           #Improved Predictors package
library(rpart.plot)      #Tree design 
library(ROCR)            #ROC curve
library(caret)           #confusion matrix: caret_6.0-76 Version
library(randomForest)    #randomForest
library(nnet)            #neural networks
library("ggplot2")
#library('RWeka')
library('ROSE')
library("optparse")
library(pROC)
library('MLmetrics')

### functions

#input: result of rfe function
filter_wilcox <-function(profile.1){
  # compare mean values of KAPPA in function of the number of variables:
  wilcox=pairwise.wilcox.test(profile.1$resample$Kappa, profile.1$resample$Variables,p.adjust.method='none') #<0.1: selection of 4 features
  ttest=pairwise.t.test(profile.1$resample$Kappa, profile.1$resample$Variables,p.adjust.method='none') #<0.1: selection of 4 features
  
  #first column with all values > 0.1 is the min number of features.
  minvalues=apply(na.omit(wilcox$p.value),2,min)
  which(minvalues>0.1)
  numbervariable=which(minvalues>0.1)[1]
  
  print("Features selected by the recursive feature elimination")
  print(profile.1$optVariables[1:numbervariable])
  return(numbervariable)
}


### definition of the inputs

option_list=list(
  make_option(c("-c","--sqanti_classif"),type="character", default = NULL, help="Sqanti classification output file"),
  make_option(c("-d","--dir"),type="character", default="Filter_out", help="Output directory name"),
  make_option(c('-w','--wilcox'),type="integer",default=0, help="Feature selection option: 0 = max rfe value, 1 = minimum number of feature with a mean non different to the mean of the classifier obtain on all the features."),
  make_option(c("-t","--pourcent_training"),type="double",default=0.8,help="the percentage of data that goes to training (parameter p of the function createDataPartition)"),
  make_option(c("-e","--exon_cov"),type="integer", default=0,help="Should the exon coverage be considered ? 0 No, 1 Yes, default 0"),
  make_option(c("-p","--TP"), type="character",default = NULL,help="file containing the list of the TP transcripts, one ID by line, no header"),
  make_option(c("-n","--TN"), type="character",default = NULL,help="file containing the list of the TN transcripts, one ID by line, no header"),
  make_option(c("-j","--threshold"), default = "0.5", help="machine learning probability threshold to classify posiive isoforms"),
  make_option(c("-i","--intrapriming"), default = "80", help="polyA percentage thereshold to flag an isoform as intra-priming")
)
 
opt_parser=OptionParser(option_list = option_list)
opt=parse_args(opt_parser) #list of the args
opt$threshold = as.numeric(opt$threshold)
opt$intrapriming = as.numeric(opt$intrapriming)

if(opt$TP=="NULL"){
  opt$TP = NULL
}

if(opt$TN=="NULL"){
  opt$TN = NULL
}
 
d <-read.table(file=opt$sqanti_classif, sep ="\t", header=TRUE,as.is =TRUE)
dme <- d[which(d$exons==1),] # store mono-exon transcript for intrapriming filtering.
d1 <- d[which(d$exons!=1),] #remove monoexon

if(nrow(d1)==0){
  print("Machine Learning filtering won't be applied because all the isoforms are mono-exon")
}else{
  print("Number of of transcripts with more than one exon:")
  print(dim(d1)[1])
}

if(is.null(opt$TP)|is.null(opt$TN)){
  print("Not provided training set. Training set from input data.")
  if(nrow(d1[d1$structural_category=="full-splice_match",])<40 | nrow(d1[which(d1$structural_category=="novel_not_in_catalog" & d1$all_canonical=="non_canonical"),])<40){
        print("Not enought Full splice match (FSM) and/or novel not in catalog-noncanonical (NNC-NC) isoforms to be used as training set. Skipping machine learning.")
        flag=FALSE
  }else{print("Full splice match as true positives, novel not in catalogue as true negatives")
        flag=TRUE}
}

d1$FSM_class=as.factor(d1$FSM_class)
d1$coding=as.factor(d1$coding)
dtotal=d1

##### ML CLASSIFIER

if(flag==TRUE){
  
  print("Data description:")
  
  #############################
  # Remove descriptive columns
  #############################
  
  if(is.null(d1$isoform)){
    rownames(d1)=d1$PacBio
    d1=d1[,-which(colnames(d1)=='PacBio')]
  }else{
    rownames(d1)=d1$isoform
    d1=d1[,-which(colnames(d1)=='isoform')]
  }
  
  colRem = c('chrom','strand','associated_gene', 'associated_transcript', 'ref_length','ref_exons','diff_to_TSS',
    'diff_to_TTS', 'ORF_length', 'CDS_length', 'CDS_start', 'CDS_end', 'perc_A_downstream_TTS')
  d3=d1[,-which(colnames(d1)%in%colRem)]
  
  
  
  #############################
  #1) Dealing with NA
  
  # if NA in Min_sample_cov, replace this value by 0
  for(i in 1:length(d3$ratio_exp)){
    if(is.na(d3$min_sample_cov[i])){
      d3$min_sample_cov[i]=0
      dtotal$min_sample_cov[i]=0
    }
  }
  
  # if NA in MinCov, replace this value by 0
  for(i in 1:length(d3$ratio_exp)){
    if(is.na(d3$min_cov[i])){
      d3$min_cov[i]=0
      dtotal$min_cov[i]=0
    }
  }
  
  # if NA in MinCovPos, replace this value by 0
  for(i in 1:length(d3$ratio_exp)){
    if(is.na(d3$min_cov_pos[i])){
      d3$min_cov_pos[i]=0
      dtotal$min_cov_pos[i]=0
    }
  }
  
  
  # if ratioExp is NA, replace this value by 0
  for(i in 1:length(d3$ratio_exp)){
    if(is.na(d3$ratio_exp[i])){
      d3$ratio_exp[i]=0
      dtotal$ratio_exp[i]=0
    }
  }
  
  # if all sdCov is NA, replace by 0. If NA is some of them, replace with median value of sdCov in the dataset
  if (all(is.na(d3$sd_cov))){
    d3$sd_cov = 2
    dtotal$sd_cov = 2
    
  }else{
    medt2=median(as.numeric(d3[!is.na(d3$sd_cov), "sd_cov"]))
    d3[is.na(d3$sd_cov),"sd_cov"] <- medt2
    dtotal[is.na(dtotal$sd_cov),"sd_cov"] <- medt2
    
  }
  
  # if NA in indels, replace this value by 0
  for(i in 1:length(d3$ratio_exp)){
    if(is.na(d3$n_indels[i])){
      d3$n_indels[i]=0
      dtotal$n_indels[i]=0
    }
  }
  
  # if NA in indels near junctions, replace this value by 0
  for(i in 1:length(d3$ratio_exp)){
    if(is.na(d3$n_indels_junc[i])){
      d3$n_indels_junc[i]=0
      dtotal$n_indels_junc[i]=0
    }
  }
  
  if(opt$exon_cov==0){
    if(!is.null(d1$FIRST_EXON_COV)){
      d3=d3[,-grep('FIRST_EXON_COV',colnames(d3))]
      d3=d3[,-grep('FIRST_EXON_NUMBER',colnames(d3))]
      d3=d3[,-grep('LAST_EXON_NUMBER',colnames(d3))]
      d3=d3[,-grep('LAST_EXON_COV',colnames(d3))]
      d3=d3[,-grep('EXON_min_cov',colnames(d3))]
    }
  }
  if(opt$exon_cov==1){
    if(is.null(d3$FIRST_EXON_COV)|is.null(d3$FIRST_EXON_NUMBER)|is.null(d3$LAST_EXON_COV)|is.null(d3$LAST_EXON_NUMBER)|is.null(d3$EXON_min_cov)){
      print("The information on the exon coverage is not available (columns missing). The variables 'FIRST_EXON_COV','FIRST_EXON_NUMBER','LAST_EXON_NUMBER,'LAST_EXON_COV','EXON_min_cov' will not be considered")
    }
  }
  
  ########################################
  # Preprocessing
  ########################################
  print("-------------------------------------------------")
  print("Preprocessing: remove of the near zero variables")
  
  nzv=nearZeroVar(d3)

  if(length(colnames(d3)[nzv])!=0){
    d27=d3[,-nzv]
    print("removed columns: ")
    print(colnames(d3)[nzv])}else{d27=d3
      print("no near zero variables! ")}
  
  print("Preprocessing: Identification of near zero variance predictors")
  dmatrix=d27
  if(!is.null(dmatrix$FSM_class)){
    dmatrix=dmatrix[,-grep('FSM_class',colnames(dmatrix))]
  }
  if(!is.null(dmatrix$structural_category)){
    dmatrix=dmatrix[,-grep('structural_category',colnames(dmatrix))]
  }
  if(!is.null(dmatrix$subcategory)){
    dmatrix=dmatrix[,-grep('subcategory',colnames(dmatrix))]
  }
  if(!is.null(dmatrix$all_canonical)){
    dmatrix=dmatrix[,-grep('all_canonical',colnames(dmatrix))]
  }
  
  
  d=data.matrix(dmatrix)
  #apply(d,2,function(x) summary(is.na(x)))
  
  descrCorr=cor(d)
  highCorr <- findCorrelation(descrCorr, cutoff=0.9)
  findCorrelation(descrCorr, cutoff=0.9,verbose = TRUE, names=TRUE)
  
  print("Preprocessing: remove of the features with a correlation > 0.9")
  
  d28=d27
  if(length(highCorr)>0){
    print("List of the features removed: ")
    for(i in highCorr){
      print(colnames(d)[i])
      d28=d27[,-grep(colnames(d)[i],colnames(d28))]
    }
  }else{
    print("No feature removed")
  }
  
  
  ###########
  # training set:
  # FSM as positives
  # NNC-NC as negatives
  ###############
  if(is.null(opt$TP)){
    print("TP: full splice match")
    print("TN: novel not in catalog non canonical")
    fsm=d28[d28$structural_category=="full-splice_match",]
    nnc=d28[d28$structural_category=="novel_not_in_catalog",]
    nncnc=nnc[nnc$all_canonical=='non_canonical',]
    trainingsetcomplet=rbind(fsm,nncnc)
    Class=factor(c(rep("POS",length(fsm$length)),rep("NEG",length(nncnc$length))))
  }else{
    TP=read.table(opt$TP,as.is=TRUE)
    TN=read.table(opt$TN,as.is=TRUE)
    Class=factor(c(rep("POS",length(TP$V1)),rep("NEG",length(TN$V1))))
    trainingset=d28[c(TP$V1,TN$V1),]
  }
  dim(trainingsetcomplet)
  trainingset=trainingsetcomplet[,-grep('structural_category',colnames(trainingsetcomplet))]
  trainingset=trainingset[,-grep('all_canonical',colnames(trainingset))]
  trainingset=trainingset[,-grep('subcategory',colnames(trainingset))]
  
  #print(summary(trainingset))
  #print(head(trainingset$coding))
  
  
  

  ####  partition
  
  print("Partition train set / test set")
  print("Proportion of the label data used for the training (between 0 and 1):")
  print(opt$pourcent_training)
  
  set.seed(3456)
  inTraining=createDataPartition(Class,p=opt$pourcent_training,list=FALSE,times=1)
  
  training=trainingset[inTraining,]
  testing=trainingset[-inTraining,]
  
  print("Description of the training set:")
  print("Table number of example positif and negatif in the training set")
  print(table(trainingsetcomplet[inTraining,]$structural_category))
  print("Table number of example positif and negatif in the test set")
  print(table(trainingsetcomplet[-inTraining,]$structural_category))
  
  
  
  ####################################
  # Machine learning
  ####################################
  
  #10 times 10 cross validation
  ctrl=trainControl(method="repeatedcv",repeats=10,
                    classProbs = TRUE,
                    summaryFunction = twoClassSummary,
                    sampling='down',returnData=TRUE,
                    savePredictions = TRUE, returnResamp='all')
  set.seed(1)
  #default: 500 trees
  #randomforest <- train(training[,profile.1$optVariables[1:numbervariable]],Class[inTraining],
  randomforest <- train(training,Class[inTraining],
                        method ='rf',
                        tuneLength=15,
                        metric = "ROC",
                        trControl = ctrl)
  
  imp=varImp(randomforest,scale=FALSE)
  write.table(imp$importance,file=paste(opt$dir,'/VarImptable.txt',sep='' ), quote=F, sep="\t")
  
  
  # Plotting Variable importance
  
  for(i in 1:length(colnames(training))){
    #distribution on the test set
    filename=paste(opt$dir,"/boxplot_",colnames(training[i]),'_training',sep='')
    # if(colnames(training)[i] %in% profile.1$optVariables[1:numbervariable]){
    #   filename=paste(filename,'_selected',sep='')
    # }
    filename=paste(filename,'.pdf',sep='')
    if(is.factor(training[,i])){
      filename=paste(opt$dir,"/hist_",colnames(training[i]),sep='')
      filename=paste(filename,'.pdf',sep='')
      pdf(filename)
      plot(training[,i]~Class[inTraining],main=colnames(training[i]),xlab="test set")
      dev.off()
    }  else{
      if(typeof(training[,i]) == "integer"|typeof(training[,i]) == "double"){
        pdf(filename)
        tmp=training[,i]+1
        boxplot(tmp~Class[inTraining],log='y',main=colnames(training[i]),xlab="test set")
        dev.off()
      }else{
        if(typeof(training[,i])=="logical"){
          filename=paste(opt$dir,"/barplot_",colnames(training[i]),sep='')
          filename=paste(filename,'.pdf',sep='')
          pdf(filename)
          plot(as.factor(training$bite)~Class[inTraining],main=colnames(training[i]),xlab="test set")
          dev.off()
        }else{
          pdf(filename)
          tmp=training[,i]+1
          boxplot(tmp~Class[-inTraining],log='y',main=colnames(training[i]),xlab="test set")
          dev.off()
        }
      }
    }
  }
  
  
  
  
  
  ###############################
  #### classifier in testing set
  ###############################
  
  test_pred_prob=predict(randomforest,testing,type='prob')
  pred=factor(ifelse(test_pred_prob$POS>=opt$threshold,"POS","NEG"))
  a=data.frame(pred,obs=Class[-inTraining],POS=test_pred_prob$POS,NEG=test_pred_prob$NEG)
  print("ROC, Sens, Spec on the test set")
  print(twoClassSummary(a,lev=levels(a$obs)))
  write("ROC, Sens, Spec on the test set",file=paste(opt$dir,'/statistics_testset.txt',sep=''),append=TRUE)
  write(twoClassSummary(a,lev=levels(a$obs)),file=paste(opt$dir,'/statistics_testset.txt',sep=''), append=TRUE)
  
  print("Area under precision-recall curve, Precision, Recall, F ")
  print(prSummary(a,lev=levels(a$obs)))
  write("Area under precision-recall curve, Precision, Recall, F ",file=paste(opt$dir,'/statistics_testset.txt',sep=''),append=TRUE)
  write(prSummary(a,lev=levels(a$obs)),file=paste(opt$dir,'/statistics_testset.txt',sep=''), append=TRUE)
  
  testpredandclass=cbind(test_pred_prob,Class[-inTraining])
  write.table(testpredandclass,paste(opt$dir,'/Pred_test_and_class.txt',sep=''))
  
  cm=confusionMatrix(data=pred,reference=Class[-inTraining],positive="POS")
  info="rows:predictions \ncolumns:reference"
  write(paste(info,"\n"),file=paste(opt$dir,"/confusion_matrix_testingset.txt",sep=''),append = TRUE)
  write.table(cm$table,file=paste(opt$dir,"/confusion_matrix_testingset.txt",sep=''),append=TRUE)
  

  ###### ROC curve
  # 1) in function of the probability on the test set (not the same proportion of positives and negatives)
  fileroc=paste(opt$dir,"/ROC_curve_ontestset.pdf",sep='')
  pdf(fileroc)
  r=roc(as.numeric(Class[-inTraining]),test_pred_prob$POS,percent = TRUE)
  auc(r)
  plot.roc(r)
  text(20,10,paste('AUC =',signif(auc(r),4)))
  text(20,5,paste('CI 95% = [',signif(ci(r)[1],4),',',signif(ci(r)[2]),']'))
  #dev.off()
  # # 2) same proportion positives and negatives on the test set:
  #testing
  #list of the testing positives
  alltestpos=which(Class[-inTraining]=='POS')
  alltestneg=which(Class[-inTraining]=='NEG')
  nbpos=length(alltestpos)
  nbneg=length(alltestneg)
  if(nbpos<nbneg){
    sampleneg=sample(alltestpos,nbpos,replace=FALSE)
    newtest=testing[c(sampleneg,alltestpos),]
    Classnewtest=factor(c(rep('NEG',nbpos),rep('POS',nbpos)))
  }else{
    samplepos=sample(alltestpos,nbneg,replace=FALSE)
    newtest=testing[c(samplepos,alltestneg),]
    Classnewtest=factor(c(rep('POS',nbneg),rep('NEG',nbneg)))
  }
  set.seed(1)
  
  test_pred_prob2=predict(randomforest,newtest,type='prob')
  
  r=roc(as.numeric(Classnewtest),test_pred_prob2$POS,percent = TRUE)
  auc(r)
  #Area under the curve: 99.29%
  #plot.roc(r)
  lines(r,col='red')
  ci(r)
  text(20,20,paste('AUC (same proportion +/-)=',signif(auc(r),4)),col='red')
  text(20,15,paste('CI 95% = [',signif(ci(r)[1],4),',',signif(ci(r)[2]),']'),col='red')
  dev.off()
  
  
  

  ###############################
  #### classifier in our dataset
  ###############################
  
  preditproba=predict(randomforest,dtotal,type='prob')
  colnames(preditproba) = gsub("NEG","NEG_MLprob", colnames(preditproba))
  colnames(preditproba) = gsub("POS","POS_MLprob", colnames(preditproba))
  
  dtotalandclassprob=cbind(dtotal,preditproba)
  
  mm.neg = dtotalandclassprob[dtotalandclassprob$POS<opt$threshold,]
  mm.neg = mm.neg[-which(mm.neg$structural_category %in% c("full-splice_match","incomplete-splice_match")),"isoform"]
  
  dtotalandclassprob[which(dtotalandclassprob$isoform %in% mm.neg),"ML_classifier"] <- "NEG" 
  dtotalandclassprob[is.na(dtotalandclassprob$ML_classifier),"ML_classifier"] <- "POS"

  
  if (nrow(dme)>0){
    dtotal = rbind(dtotalandclassprob, data.frame(dme, NEG_MLprob=NA, POS_MLprob=NA, ML_classifier="NA"))
    print("ML classifier results (NEG: Isoforms classified as Artifacts)")
    print(table(dtotal$ML_classifier))
  }
} else{dtotal=data.frame(d, ML_classifier="NA")}



##### INTRA-PRIMING FILTERING

dtotal[,"intra-priming"] = dtotal$perc_A_downstream_TTS > as.numeric(opt$intrapriming) & !(dtotal$structural_category %in% c("full-splice_match","incomplete-splice_match") )

print("Intra-priming filtered transcripts:")
print(table(dtotal$`intra-priming`))
  


##### WRITING RESULTS

dtotal[which(dtotal$ML_classifier=="NEG" | dtotal$`intra-priming`==TRUE),"SQANTI_filter"] <- "Artifact"
dtotal[is.na(dtotal$SQANTI_filter),"SQANTI_filter"] <- "Isoform"
fileres=paste(opt$dir,'/',basename(opt$sqanti_classif),"_filterResults.txt",sep='')
write.table(dtotal,fileres,sep='\t',quote = F, row.names = F)

fileres2=paste(opt$dir,'/',basename(opt$sqanti_classif),"_curatedTranscriptome.txt",sep='')
write.table(dtotal[which(dtotal$SQANTI_filter=="Isoform"),"isoform"],fileres2,sep='\t',quote = F, row.names = F, col.names = F)

print("*** SQANTI filter results:")
print(table(dtotal$SQANTI_filter))

print("*** SQANTI filtering finished!")






