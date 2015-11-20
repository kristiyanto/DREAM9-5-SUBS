# THIS FILE LOADS AND CLEANS THE CORETABLE_TEST_KY.CSV FILE AND SAVES A SUBMISSION.CSV FILE THAT CONTAINS
# THE RPT VALUE, THE RISK, AND THE PREDICTION OF THAT PATIENT USING THE MODEL SAVED IN AZURE ML
## Danny only: 
#setwd("~/Dropbox/588/DREAM PROJECT/FINAL RESULTS/AFTER FINAL")


require(FSelector) 
require(e1071)
require(randomForest)
require(timeROC)
library(FSelector) 
library(randomForest)
library(e1071)
library(impute)
library(nnet)

source("inc_datacleanup.r")
source("score.R")

### READ IN DATA AND SPLIT THE DATA
table.for.model   <- read.csv("CoreTable_training.csv", stringsAsFactors = F)
table.for.testing <- read.csv("clean_validation.csv", stringsAsFactors = F)
table.for.testing <- cleanData(table.for.testing)
split.table       <- split(table.for.model, table.for.model$STUDYID)
summary           <- NULL
evaluation        <- NULL
contingency.table <- NULL
selected.features <- NULL

## 3 fold validation
for (fold in c(1,2,3))  ## LOOP THROUGH THE DATA FRAMES
{
  ## PICK THE TRAINING AND TESTING COLUMNS DEPENDING ON THE FOLD.
  if (fold == 1)
  {
    training.cols = c(2,3)
    testing.col   = c(1) 
    testing.name  = "ACCENT2"
  } else if (fold == 2) {
    training.cols  = c(1,3)
    testing.col    = c(2) 
    testing.name   = "CELGENE"
  } else {
    training.cols = c(1,2)
    testing.col   = c(3) 
    testing.name  = "EFC6546"
  }
  

  curr.training.data <- rbind(split.table[[training.cols[1]]], split.table[[training.cols[2]]])
  curr.testing.data  <- split.table[[testing.col]]
  
  
  
  curr.training.data <- cleanData(curr.training.data)
  curr.testing.data  <- cleanData(curr.testing.data)
  
  ## Double Check for missing values
  #sapply(curr.training.data, function(x) sum(is.na(x)))
  
  #### METADATA EDITOR
  all.features          <- intersect(colnames(curr.training.data),colnames(curr.testing.data))
  factorcols <- c('STUDYID',	'RPT',	'DISCONT',	'GLEAS_DX',	'TSTAG_DX',	'RACE_C',	'REGION_C',	'SMOKE',	'SMOKFREQ',	
                  'SMOKSTAT',	'ECOG_C',	'TRT2_ID',	'TRT3_ID',	'NON_TARGET',	'TARGET',	'BONE',	'RECTAL',	'LYMPH_NODES',	
                  'KIDNEYS',	'LUNGS',	'LIVER',	'PLEURA',	'OTHER',	'PROSTATE',	'ADRENAL',	'BLADDER',	'PERITONEUM',	
                  'COLON',	'SOFT_TISSUE',	'ABDOMINAL',	'ORCHIDECTOMY',	'PROSTATECTOMY',	'TURP',	'LYMPHADENECTOMY',	
                  'SPINAL_CORD_SURGERY',	'BILATERAL_ORCHIDECTOMY',	'PRIOR_RADIOTHERAPY',	'ANALGESICS',	'ANTI_ANDROGENS',	
                  'GLUCOCORTICOID',	'GONADOTROPIN',	'BISPHOSPHONATE',	'CORTICOSTEROID',	'IMIDAZOLE',	'ACE_INHIBITORS',	
                  'BETA_BLOCKING',	'HMG_COA_REDUCT',	'ESTROGENS',	'ANTI_ESTROGENS',	'ARTTHROM',	'CEREBACC',	'CHF',	
                  'DVT',	'DIAB',	'GASTREFL',	'GIBLEED',	'MI',	'PUD',	'PULMEMB',	'PATHFRAC',	'SPINCOMP',	'COPD',	
                  'MHBLOOD',	'MHCARD',	'MHCONGEN',	'MHEAR',	'MHENDO',	'MHEYE',	'MHGASTRO',	'MHGEN',	'MHHEPATO',	
                  'MHIMMUNE',	'MHINFECT',	'MHINJURY',	'MHINVEST',	'MHMETAB',	'MHMUSCLE',	'MHNEOPLA',	'MHNERV',	'MHPSYCH',	
                  'MHRENAL',	'MHRESP',	'MHSKIN',	'MHSOCIAL',	'MHSURG',	'MHVASC')
  all.features <- intersect(all.features, factorcols)
  for(i in 1:length(all.features))
  {
    cname <- all.features[i]
    level.list <- union(levels(curr.training.data[,cname]), levels(curr.testing.data[,cname]))
    levels(curr.training.data[,cname])  <- level.list
    levels(curr.testing.data[,cname])   <- level.list
  }

  ## BAGGING 
  ### Split the data into 2 based on DISCONT
  bagA      <- subset(curr.training.data, DISCONT==0)
  bagB      <- subset(curr.training.data, DISCONT==1)
  set.seed(2)
  #bagA      <- bagA[sample(nrow(bagA), round(0.5*nrow(bagA)), replace = F),]
  bagA      <- bagA[sample(nrow(bagA), round(0.5*nrow(bagA)), replace = F),]
  curr.training.data <- as.data.frame(rbind(bagA,bagB))
  plot(curr.training.data$DISCONT,main=testing.name)
  

  ## FEATURE SELECTION
  ##### Run ChiSquared
  col_names              <- names(curr.training.data) 
  table.descrete         <- lapply(curr.training.data[,col_names] , factor)
  table.descrete$RPT     <- NULL
  table.descrete$STUDYID <- NULL
  weight                 <- chi.squared(DISCONT ~ ., data=table.descrete)

  subsetChi             <- c("DISCONT", cutoff.k(weight,50))
  
  #### Random Forest Importance
  table.descrete        <- curr.training.data
  table.descrete$RPT    <- NULL
  table.descrete$STUDYID        <- NULL
  weight                <- random.forest.importance(DISCONT ~ ., data=table.descrete)
  subsetRF              <- c("DISCONT", cutoff.k.percent(weight,80))
  

  #### Removing Features that are not in the Testing Data
  subsetRF              <- intersect(subsetRF, colnames(curr.testing.data))
  subsetChi             <- intersect(subsetChi, colnames(curr.testing.data))
  
  evaluatornb <- function(subsetHill) {
    k <- 1
    results = sapply(1:k, function(i) {
      subsetHill <- c(subsetHill, "DISCONT")
      model.nb        <- naiveBayes(DISCONT ~., data=curr.training.data[,subsetHill],probability=T)
      result.nb        <- predict(model.nb,   curr.testing.data, type = "raw")
      probs.nb        <- result.nb[ ,2]
      nb.auc <- score_q2(probs.nb , curr.testing.data$DISCONT)
      return(nb.auc)
    })
    print(subsetHill)
    print(mean(results))
    return(mean(results))
  }
  
  # TEMP
  #subsetHill <- hill.climbing.search(colnames((curr.training.data)[-c(1,2,3)]), evaluatornb)
  
  # Instead of rerun the brute feature selection, load from file
  subsetHill <- read.csv(file = paste("SelectedFeatures",testing.name,".csv"))
  subsetHill <- subsetHill$x
  
  subsetHill <- intersect(subsetHill,colnames(curr.testing.data))

  #write.csv(subsetHill, file = paste("SelectedFeatures",testing.name,".csv"))
  ## BUILD THE MODELS (TRAINING)
  
  model.rf        <- randomForest(DISCONT ~., data=curr.training.data[,subsetHill], na.action = na.omit, ntree=28, probability=T)
  model.svm       <- svm(DISCONT ~., data=curr.training.data[,subsetChi], kernel = "linear", probability=T )
  model.nnet      <- nnet(DISCONT ~., data=curr.training.data[,subsetChi], size=50, rang=0.1,maxit=160,MaxNWts=5000)
  model.nb        <- naiveBayes(DISCONT ~., data=curr.training.data[,subsetHill],probability=T)
  
  
  ## RUN THE PREDICTION

  result.rf       <- predict(model.rf,   curr.testing.data, type = "prob")
  result.svm      <- predict(model.svm,  curr.testing.data, probability =T)
  result.nnet     <- predict(model.nnet, curr.testing.data, type = "raw")
  result.nb       <- predict(model.nb,   curr.testing.data, type = "raw")


  
  
  ## PUT THE RESULT TOGETHER (VOTING BUCKET)
  probs.rf        <- result.rf[ ,2]
  probs.svm       <- attr(result.svm, "probabilities")[ ,2]
  probs.nb        <- result.nb[ ,2]
  
  
  result.ens      <- (probs.rf   * 1 +
                      probs.svm  * 1 +
                      probs.nb   * 1 +
                      result.nnet * 1) / 3
  
  foldsummary    <- cbind(fold,curr.testing.data[,c("RPT","STUDYID")],probs.rf, probs.svm, result.nnet, probs.nb, result.ens, curr.testing.data$DISCONT)

  #taking out NB for a moment
 # foldsummary    <- cbind(fold,curr.testing.data[,c("RPT","STUDYID")],probs.rf, probs.svm, result.nnet, NA, result.ens, curr.testing.data$DISCONT)
  
   colnames(foldsummary) <- c("Fold","RPT","TestingStudy","Probs.RF","Probs.SVM", "Probs.NNET","Probs.NB", "Results.Ensemble", "ACTUAL")
   summary        <- rbind(summary,foldsummary)

   
   ###  USE THE SUMMARY TO SCORE THE MODEL
   
   #RF.AUC
   rf.auc <- score_q2(foldsummary$Probs.RF  , foldsummary$ACTUAL)
   
   #RF.ACC | #RF.F1
   TN <- sum(round(foldsummary$Probs.RF) == 0 & foldsummary$ACTUAL == 0)
   TP <- sum(round(foldsummary$Probs.RF) == 1 & foldsummary$ACTUAL == 1)  
   FP <- sum(round(foldsummary$Probs.RF) == 1 & foldsummary$ACTUAL == 0)  
   FN <- sum(round(foldsummary$Probs.RF) == 0 & foldsummary$ACTUAL == 1)
   prec <- TP / (TP + FP)
   rec  <- TP / (TP + FN)
   
   rf.F1  <- (2 * prec * rec) / (prec + rec)
   rf.acc <- (TN   + TP)  / (TN + TP + FP + FN)
   
   contingency.table <- rbind(contingency.table, c(fold, "RF", TP, TN, FP, FN))
   
   #SVM.AUC
   svm.auc <- score_q2(foldsummary$Probs.SVM  , foldsummary$ACTUAL)
   #SVM.ACC | #SVM.F1
   TN <- sum(round(foldsummary$Probs.SVM) == 0 & foldsummary$ACTUAL == 0)
   TP <- sum(round(foldsummary$Probs.SVM) == 1 & foldsummary$ACTUAL == 1)  
   FP <- sum(round(foldsummary$Probs.SVM) == 1 & foldsummary$ACTUAL == 0)  
   FN <- sum(round(foldsummary$Probs.SVM) == 0 & foldsummary$ACTUAL == 1)
   prec <- TP / (TP + FP)
   rec  <- TP / (TP + FN)
   
   svm.F1  <- (2 * prec * rec) / (prec + rec)
   svm.acc <- (TN   + TP)  / (TN + TP + FP + FN)
   contingency.table <- rbind(contingency.table, c(fold, "SVM", TP, TN, FP, FN))
   
   #NNET.AUC
   nnet.auc <- score_q2(foldsummary$Probs.NNET  , foldsummary$ACTUAL)
   
   #NNET.ACC | #NNET.F1
   TN <- sum(round(foldsummary$Probs.NNET) == 0 & foldsummary$ACTUAL == 0)
   TP <- sum(round(foldsummary$Probs.NNET) == 1 & foldsummary$ACTUAL == 1)  
   FP <- sum(round(foldsummary$Probs.NNET) == 1 & foldsummary$ACTUAL == 0)  
   FN <- sum(round(foldsummary$Probs.NNET) == 0 & foldsummary$ACTUAL == 1)
   prec <- TP / (TP + FP)
   rec  <- TP / (TP + FN)
   
   nnet.F1  <- (2 * prec * rec) / (prec + rec)
   nnet.acc <- (TN   + TP)  / (TN + TP + FP + FN)
   contingency.table <- rbind(contingency.table, c(fold, "NNET", TP, TN, FP, FN))
   
   
   #NB.AUC
   nb.auc <- score_q2(foldsummary$Probs.NB  , foldsummary$ACTUAL)
   #NB.ACC | #NB.F1
   TN <- sum(round(foldsummary$Probs.NB) == 0 & foldsummary$ACTUAL == 0)
   TP <- sum(round(foldsummary$Probs.NB) == 1 & foldsummary$ACTUAL == 1)  
   FP <- sum(round(foldsummary$Probs.NB) == 1 & foldsummary$ACTUAL == 0)  
   FN <- sum(round(foldsummary$Probs.NB) == 0 & foldsummary$ACTUAL == 1)
   prec <- TP / (TP + FP)
   rec  <- TP / (TP + FN)
   
   nb.F1  <- (2 * prec * rec) / (prec + rec)
   nb.acc <- (TN   + TP)  / (TN + TP + FP + FN)
   contingency.table <- rbind(contingency.table, c(fold, "NB", TP, TN, FP, FN))
   
   #ENS.AUC
   ens.auc <- score_q2(foldsummary$Results.Ensemble  , foldsummary$ACTUAL)
   #NB.ACC | #NB.F1
   TN <- sum(round(foldsummary$Results.Ensemble) == 0 & foldsummary$ACTUAL == 0)
   TP <- sum(round(foldsummary$Results.Ensemble) == 1 & foldsummary$ACTUAL == 1)  
   FP <- sum(round(foldsummary$Results.Ensemble) == 1 & foldsummary$ACTUAL == 0)  
   FN <- sum(round(foldsummary$Results.Ensemble) == 0 & foldsummary$ACTUAL == 1)  
   prec <- TP / (TP + FP)
   rec  <- TP / (TP + FN)
   
   ens.F1  <- (2 * prec * rec) / (prec + rec)
   ens.acc <- (TN   + TP)  / (TN + TP + FP + FN)
   contingency.table <- rbind(contingency.table, c(fold, "ENS", TP, TN, FP, FN))
   
# pulling out NB for  amoment
  evaluation.row <- cbind(fold, rf.auc, svm.auc, nnet.auc, nb.auc, ens.auc
                          , rf.acc, svm.acc, nnet.acc, nb.acc, ens.acc
                          , rf.F1,  svm.F1, nnet.F1, nb.F1,  ens.F1)
   
#    evaluation.row <- cbind(fold, rf.auc, svm.auc, nnet.auc, NA, ens.auc
#                            , rf.acc, svm.acc, nnet.acc, NA, ens.acc
#                            , rf.F1,  svm.F1, nnet.F1, NA,  ens.F1)
#    
   
   colnames(evaluation.row) <- c("Fold","AUC.RF", "AUC.SVM", "AUC.NNET", "AUC.NB", "AUC.Ens"
                             , "ACC.RF", "ACC.SVM", "ACC.NNET", "ACC.NB", "ACC.Ens"
                             , "F1.RF" , "F1.SVM" , "F1.NNET", "F1.NB" , "F1.Ens")  
   evaluation <- rbind(evaluation, evaluation.row)
} 

rownames(evaluation) <- c("ACCENT2","CELGENE", "EFC6546")


#REMOVING NB SO THERE ARE ONLY FOUR MODELS INSTEAD OF 5
rownames(contingency.table) <- c(rep("ACCENT2",5), rep("CELGENE",5), rep("EFC6546",5))
# rownames(contingency.table) <- c(rep("ACCENT2",4), rep("CELGENE",4), rep("EFC6546",4))
colnames(contingency.table) <- c("Fold", "Model", "TP","TN","FP","FN")


evaluation
write.csv(evaluation, file="Evaluation.csv")
write.csv(contingency.table, file = "Contingency.csv")

result <- predict(model.nb, table.for.testing, type="probs")
