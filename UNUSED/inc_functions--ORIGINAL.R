
#### PREDEDINED COLUMN TYPES
factor.cols   <- c('STUDYID',	'GLEAS_DX',	'TSTAG_DX',	'RACE_C',	'REGION_C',	'SMOKE',	'SMOKFREQ',	
                    'SMOKSTAT',	'ECOG_C',	'TRT2_ID',	'TRT3_ID',	'NON_TARGET',	'TARGET',	'BONE',	'RECTAL',	'LYMPH_NODES',	
                    'KIDNEYS',	'LUNGS',	'LIVER',	'PLEURA',	'OTHER',	'PROSTATE',	'ADRENAL',	'BLADDER',	'PERITONEUM',	
                    'COLON',	'SOFT_TISSUE',	'ABDOMINAL',	'ORCHIDECTOMY',	'PROSTATECTOMY',	'TURP',	'LYMPHADENECTOMY',	
                    'SPINAL_CORD_SURGERY',	'BILATERAL_ORCHIDECTOMY',	'PRIOR_RADIOTHERAPY',	'ANALGESICS',	'ANTI_ANDROGENS',	
                    'GLUCOCORTICOID',	'GONADOTROPIN',	'BISPHOSPHONATE',	'CORTICOSTEROID',	'IMIDAZOLE',	'ACE_INHIBITORS',	
                    'BETA_BLOCKING',	'HMG_COA_REDUCT',	'ESTROGENS',	'ANTI_ESTROGENS',	'ARTTHROM',	'CEREBACC',	'CHF',	
                    'DVT',	'DIAB',	'GASTREFL',	'GIBLEED',	'MI',	'PUD',	'PULMEMB',	'PATHFRAC',	'SPINCOMP',	'COPD',	
                    'MHBLOOD',	'MHCARD',	'MHCONGEN',	'MHEAR',	'MHENDO',	'MHEYE',	'MHGASTRO',	'MHGEN',	'MHHEPATO',	
                    'MHIMMUNE',	'MHINFECT',	'MHINJURY',	'MHINVEST',	'MHMETAB',	'MHMUSCLE',	'MHNEOPLA',	'MHNERV',	'MHPSYCH',	
                    'MHRENAL',	'MHRESP',	'MHSKIN',	'MHSOCIAL',	'MHSURG',	'MHVASC', 'CancerSums', 'MedBodySums')

binary.cols    <- c('NON_TARGET',	'TARGET',	'BONE',	'RECTAL',	'LYMPH_NODES',	'KIDNEYS',	'LUNGS',	'LIVER',	'PLEURA',	
                 'OTHER',	'PROSTATE',	'ADRENAL',	'BLADDER',	'PERITONEUM',	'COLON',	'HEAD_AND_NECK',	'SOFT_TISSUE',	
                 'STOMACH',	'PANCREAS',	'THYROID',	'ABDOMINAL',	'ORCHIDECTOMY',	'PROSTATECTOMY',	'TURP',	'LYMPHADENECTOMY',	
                 'SPINAL_CORD_SURGERY',	'BILATERAL_ORCHIDECTOMY',	'PRIOR_RADIOTHERAPY',	'ANALGESICS',	'ANTI_ANDROGENS',	
                 'GLUCOCORTICOID',	'GONADOTROPIN',	'BISPHOSPHONATE',	'CORTICOSTEROID',	'IMIDAZOLE',	'ACE_INHIBITORS',	
                 'BETA_BLOCKING',	'HMG_COA_REDUCT',	'ESTROGENS',	'ANTI_ESTROGENS',	'ARTTHROM',	'CEREBACC',	'CHF',	'DVT',	
                 'DIAB',	'GASTREFL',	'GIBLEED',	'MI',	'PUD',	'PULMEMB',	'PATHFRAC',	'SPINCOMP',	'COPD',	'MHBLOOD',	'MHCARD',	
                 'MHCONGEN',	'MHEAR',	'MHENDO',	'MHEYE',	'MHGASTRO',	'MHGEN',	'MHHEPATO',	'MHIMMUNE',	'MHINFECT',	'MHINJURY',	
                 'MHINVEST',	'MHMETAB',	'MHMUSCLE',	'MHNEOPLA',	'MHNERV',	'MHPSYCH',	'MHRENAL',	'MHRESP',	'MHSKIN',	'MHSOCIAL',	
                 'MHSURG',	'MHVASC','ClustLab', 'ClustMedHistory')
numeric.cols  <- c('ALP',	'ALT',	'AST',	'CA',	'CREAT',	'HB',	'LDH',	'NEU',	'PLT',	'PSA',	'TBILI',	'TESTO',	'WBC',	
                  'CREACL',	'NA',	'MG',	'PHOS',	'ALB',	'TPRO',	'RBC',	'LYM',	'BUN',	'CCRC',	'GLU',	'CREACLCA')


## FUNCTION AS USED FOR GREEDY METHOD IN FEATURE SELECTIONS

evaluatornb <- function(subset.feature) {
  k <- 1                                            # Cross-fold validation
  results = sapply(1:k, function(i) {
    bag0.subset   <- bag0[sample(nrow(bag0), 
                                 round(nrow(bag0)*baggingRatio), replace = F),]
    tmp             <- as.data.frame(rbind(bag1, bag0.subset))
    model.nb        <- naiveBayes(DISCONT ~., data=tmp[,c("DISCONT",subset.feature)],probability=T)
    result.nb       <- predict(model.nb, curr.testing.data, type = "raw")[,2]
    nb.auc <- score_q2(result.nb , curr.testing.data$DISCONT)
    return(nb.auc)
    tmp             <- NULL
  })
  print(paste("NB, Mean (",k,"cv) AUC:",mean(results)))
  return(mean(results))
}

evaluatorrf <- function(subset.feature) {
  k <- 3
  results = sapply(1:k, function(i) {
    bag0.subset     <- bag0[sample(nrow(bag0), 
                                   round(nrow(bag0)*baggingRatio), replace = F),]
    tmp             <- as.data.frame(rbind(bag1,bag0.subset))
    model.rf        <- randomForest(DISCONT ~., data=tmp[,c("DISCONT",subset.feature)], mtry=round(length(subset.feature)/3), na.action = na.omit, probability=T)
    result.rf       <- predict(model.rf,   curr.testing.data, type = "prob")[,2]
    rf.auc <- score_q2(result.rf , curr.testing.data$DISCONT)
    return(rf.auc)
  })
  print(paste("RF, FS - Climbing with current AUC:",mean(results)))
  return(mean(results))
}

evaluatorsvm <- function(subset.feature) {
  k <- 3
  results = sapply(1:k, function(i) {
    bag0.subset     <- bag0[sample(nrow(bag0), 
                                   round(nrow(bag0)*baggingRatio), replace = F),]
    tmp             <- as.data.frame(rbind(bag1,bag0.subset))
    model.svm       <- svm(DISCONT ~., data=tmp[,c("DISCONT",subset.feature)], kernel="linear", scale = T)
    result.svm      <- as.factor(predict(model.svm,   curr.testing.data))
    svm.acc <- (sum(curr.testing.data$DISCONT==result.svm))/nrow(curr.testing.data)
    #rf.auc <- score_q2(probs.rf , tmp2$DISCONT)
    return(svm.acc)
  })
  print(paste("SVM",mean(results)))
  return(mean(results))
}


evaluatornnet <- function(subset.feature) {
  k <- 3
  results = sapply(1:k, function(i) {
    bag0.subset     <- bag0[sample(nrow(bag0), 
                                   round(nrow(bag0)*baggingRatio), replace = F),]
    tmp             <- as.data.frame(rbind(bag1,bag0.subset))
    model.nnet      <- nnet(DISCONT ~., data=tmp[,c("DISCONT",subset.feature)], size=length(subset.feature), trace=F)
    result.nnet       <- predict(model.nnet,   curr.testing.data, type = "raw")
    #nnet.auc <- score_q2(result.nnet , curr.testing.data$DISCONT)
    nnet.auc <- (sum(curr.testing.data$DISCONT==(round(result.nnet))))/nrow(curr.testing.data)
    return(nnet.auc)
  })
  print(paste("nnet",mean(results)))
  return(mean(results))
}

get.training.data <- function(fold)
{
  if (fold == 1)
  {
    training.cols = c(4,3)
    testing.col   = c(1) 
    testing.name  = "ACCENT2"
    
  } else if (fold == 2) {
    training.cols  = c(1,3)
    testing.col    = c(4) 
    testing.name   = "EFC6546"
    
  } else {
    training.cols = c(1,4)
    testing.col   = c(3) 
    testing.name  = "CELGENE"
  }
  
  ## Split between training and testing
  curr.training.data <- as.data.frame(rbind(split.table[[training.cols[1]]], split.table[[training.cols[2]]]))
  curr.testing.data  <- as.data.frame(split.table[[testing.col]])
  return(curr.training.data)
}

get.testing.data <- function(fold)
{
  if (fold == 1)
  {
    training.cols = c(4,3)
    testing.col   = c(1) 
    testing.name  = "ACCENT2"
    
  } else if (fold == 2) {
    training.cols  = c(1,3)
    testing.col    = c(4) 
    testing.name   = "EFC6546"
    
  } else {
    training.cols = c(1,4)
    testing.col   = c(3) 
    testing.name  = "CELGENE"
  }
  
  ## Split between training and testing
  curr.training.data <- as.data.frame(rbind(split.table[[training.cols[1]]], split.table[[training.cols[2]]]))
  curr.testing.data  <- as.data.frame(split.table[[testing.col]])
  return(curr.testing.data)
}

get.testing.name <- function(fold)
{
  if (fold == 1)
  {
    testing.name  = "ACCENT2"
  } else if (fold == 2) {
    testing.name   = "EFC6546"
    
  } else {
    testing.name  = "CELGENE"
  }
  return(testing.name)
}