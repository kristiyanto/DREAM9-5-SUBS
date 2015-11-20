
require(FSelector) 
require(e1071)
require(timeROC)
require(caret)
require(unbalanced)
library(FSelector) 
library(e1071)
library(impute)
library(caret)
library(randomForest)


source("inc_datacleanup.r")
source("score.R")
source("inc_functions.R")

## GLOBAL VARIABLES AND PARAMETERS
baggingRatio          = 0.9    # Used by Feature Selection Evaluation.
crossfold             = 10     # Number of Bootstrap to tune model parameter
FOLD                  = c(1,2,3) # Itterate through every study with following target 
                                 # (and use the other two as training data) 
                                 # 1. ACCENT2, 2. EFC6546, 3. CELGENE 

## GET THE DATA
core_table_training   <- read.csv("CoreTable_training.csv"  , stringsAsFactors = F)
core_table_validation <- read.csv("CoreTable_validation.csv", stringsAsFactors = F)

## GET THE DATA CLEANED UP
#core_table_validation$DISCONT <- 0
combined_data                 <- rbind(core_table_training, core_table_validation)
clean_data                    <- dream9.cleanData(combined_data)

## SOME DATA CLEAN UP
clean_data$DEATH[clean_data$DEATH == "."]           <- NA
clean_data$DEATH[clean_data$DEATH == "YES"]         <- 1
clean_data$DEATH[clean_data$DEATH == ""]            <- 0

clean_data$LKADT_P[clean_data$LKADT_P == "."]       <- NA
clean_data$ENDTRS_C[clean_data$ENDTRS_C == "."]     <- "UNKNOWN"
clean_data$ENTRT_PC[clean_data$ENTRT_PC == "."]     <- NA
clean_data$PER_REF[clean_data$PER_REF == "."]       <- NA
clean_data$LKADT_REF[clean_data$LKADT_REF == "."]   <- NA
clean_data$ENTRT_PC[clean_data$ENTRT_PC == "."]     <- NA


## SOME METADATA EDITOR
dependant.cols                <- c("DEATH", "LKADT_P", "ENDTRS_C", "ENTRT_PC", "PER_REF", "LKADT_REF","GLEAS_DX","TSTAG_DX")
clean_data                    <- clean_data[,!(colnames(clean_data) %in% dependant.cols)]

## VARIABLE EXTRAPOLATE
to.explode                    <- c("AGEGRP2","RACE_C","REGION_C","SMOKE","SMOKFREQ","SMOKSTAT","ECOG_C","TRT2_ID","TRT3_ID")
tmp                           <- clean_data[,to.explode]
f                             <- paste("~0+",paste(to.explode, collapse = "+"))
tmp                           <- model.matrix(as.formula(f), data=tmp)
tmp                           <- apply(tmp,2,factor)
clean_data                    <- cbind(clean_data,tmp)
to.drop                       <- c(to.explode)
clean_data                    <- clean_data[,!(colnames(clean_data) %in% to.drop)]
row.names(clean_data)         <- clean_data$RPT

## SPLIT THE DATA BACK TO TRAINING AND VALIDATION
table.for.validation          <- clean_data[clean_data$STUDYID == "AZ", ]
table.for.model               <- clean_data[((clean_data$STUDYID != "AZ") & !is.na(clean_data$DISCONT)), ]

## DROP DEPENDANT VARIABLES AWAY
to.drop                     <- c("LKADT_P","DEATH","PER_REF","LKADT_REF","LKADT_PER","GLEAS_DX", "TSTAG_DX")
table.for.model             <- table.for.model[,!(colnames(table.for.model) %in% to.drop)]

table.for.validation        <- clean_validation
split.table                 <- split(table.for.model, table.for.model$STUDYID)

## SOME VARIABLES
SCORING.TABLE               <- NULL
OUTPUT.TABLE                <- NULL


for (fold in FOLD)  ## LOOP THROUGH THE DATA FRAMES
{
  ## RETRIEVE THE CLEAN DATA SET 
  ## THE FUNCTION ASSUMES THAT THE CORE DATA IS ALREADY LOADED
  curr.training.data <- as.data.frame(get.training.data(fold=fold))
  curr.testing.data  <- as.data.frame(get.testing.data(fold=fold))
  testing.name       <- get.testing.name(fold)
  
  ## DROP VARIABLES THAT ARE NOT AVAILABLE IN TESTING FROM TRAINING AND SOME REDUDANT
  to.drop <- union(to.drop, c("SMOKFREQ",	"SMOKSTAT", "HEIGHTBL", "WEIGHTBL",	"WEIGHT", "NON_TARGET", "AGE", "X"))
  tmp     <- curr.testing.data[curr.testing.data$DISCONT==1,]
  tmp     <- tmp[,-c(1,2,3)]
  tmp     <- sapply(tmp, var)
  tmp     <- tmp[tmp==0]
  to.drop <- union(to.drop, attr(tmp,"names"))
  curr.training.data <- curr.training.data[,!names(curr.training.data) %in% to.drop]
  
  
  ## UPDATE THE LIST OF VARIABLE GROUPS
  all.features          <- intersect(colnames(curr.training.data),colnames(curr.testing.data))
  binary.cols           <- intersect(all.features, binary.cols)
  numeric.cols          <- intersect(all.features, numeric.cols)
  cols.to.convert       <- intersect(all.features, factor.cols)
  
  for(i in 1:length(binary.cols))
  {
    cname                       <- binary.cols[i]
    curr.training.data[,cname]  <- as.numeric(curr.training.data[,cname])
  }
  
  ## BALANCE THE TRAINING DATA
  library(mlr)
  library(unbalanced)
  balancer.ubunder                <- ubUnder(X=curr.training.data[,-c(3)], Y=curr.training.data[,'DISCONT'], perc = 40, method = "percUnder")
  curr.training.data              <- cbind(balancer.ubunder$Y, balancer.ubunder$X)
  names(curr.training.data)[1]    <- paste("DISCONT")
  detach(package:unbalanced, unload=T, force=T) 
  detach(package:mlr, unload=T, force=T)
  
  
  ### METADATA EDITOR
  for(i in 1:length(cols.to.convert))
  {
    cname                               <- cols.to.convert[i]
    # Convert to R Valid Names. Caret is very picky
    curr.testing.data[,cname]           <- as.factor(make.names(curr.testing.data[,cname])) 
    curr.training.data[,cname]          <- as.factor(make.names(curr.training.data[,cname]))
    # Make sure that factor levels between training and testing match
    level.list <- union(levels(curr.training.data[,cname]), levels(curr.testing.data[,cname]))
    levels(curr.training.data[,cname])  <- level.list
    levels(curr.testing.data[,cname])   <- level.list
  }
  
  
  #### FEATURE SELECTION 
  ## A SUBSET REQUIRED BY FEATURE SELECTION EVALUATION
  bag0          <- subset(curr.training.data, DISCONT=="0")
  bag1          <- subset(curr.training.data, DISCONT=="1")
  
  set.seed(123)
  features              <- setdiff(all.features,union(numeric.cols,binary.cols))
  features              <- features[-c(1,2,3)]
  subset.main           <- hill.climbing.search(features, evaluatorrf)
  subset.main           <- union(subset.main,numeric.cols)
  subset.main           <- hill.climbing.search(subset.main, evaluatorrf)
  
  print("Hill Climbing with binary.")
  subset.main           <- union(subset.main,binary.cols)
  subset.main           <- hill.climbing.search(subset.main, evaluatorrf)
  features              <- subset.main
  
  library(caret)
  set.seed(123)
  caret.training.data   <- (curr.training.data[,features])
  target                <- as.factor(make.names(curr.training.data$DISCONT,unique = F))   
  
  model.control         <- trainControl(method = "cv", number=crossfold)
  grid.rf               <- expand.grid(mtry = seq(round(ncol(caret.training.data)/10),ncol(caret.training.data),round(ncol(caret.training.data)/10)))
  print("Training is in session. Please be patient.")

  model.rf              <- train(x=caret.training.data, y=target, method = 'rf', trControl = model.control, tuneGrid = grid.rf)
  
  
  prob <- predict(model.rf, curr.testing.data, type = "prob")[,2]
  
  table.prob <- cbind(curr.testing.data[,features],prob, as.factor(round(prob)), curr.testing.data$DISCONT)
  SCORING.TABLE <- rbind(SCORING.TABLE, c(testing.name, dream9.score(prob, curr.testing.data$DISCONT)))
  ACC = round(score_q2(prob, curr.testing.data$DISCONT)*100)
  write.csv(table.prob, file = paste(testing.name,"-",ACC,".csv", sep=""))
}

print(SCORING.TABLE)

