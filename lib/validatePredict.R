# validatePredict
# Description: This library contains the function utilized to generate perform a k-fold validation 
# and subsequent calculation of prediction metrics (RMSE, MAPE, MAD) for PLS & LM

#Function for Generating and Evaluating Out-of-sample predictions using 10-fold cross-validation
validatePredict <- function(dataset, smMatrix, mmMatrix, maxIt=300, stopCriterion=7,kfold=10, reps = 10){
  
  #Identify variables to be tested
  uniqueTarget <- unique(smMatrix[,2])
  items <- NULL
  for (i in 1:length(uniqueTarget)){
    items <- c(items, mmMatrix[mmMatrix[, "latent"] == uniqueTarget[i],"measurement"])
  }
  uniqueSource <- unique(smMatrix[,1])
  sources <- NULL
  for (i in 1:length(uniqueSource)){
    sources <- c(sources,mmMatrix[mmMatrix[,"latent"]==uniqueSource[i],"measurement"])
  }
  lmtarget <- ifelse(length(intersect(uniqueTarget, uniqueSource)) == 0, uniqueTarget,setdiff(uniqueTarget, uniqueSource))  
  targets <- NULL
  for (i in 1:length(lmtarget)){
    targets <- c(targets, mmMatrix[mmMatrix[, "latent"] == lmtarget[i],"measurement"])
  }
  
  # Initialize matrices for prediction metrics
  # Initialize RMSE holders
  PLSRMSE <- matrix(,nrow=reps,ncol=length(targets),byrow =TRUE,dimnames = list(1:reps,targets))
  PLSSSE <- matrix(,nrow=kfold,ncol=length(targets),byrow =TRUE,dimnames = list(1:kfold,targets))
  LMRMSE <- matrix(,nrow=reps,ncol=length(targets),byrow =TRUE,dimnames = list(1:reps,targets))
  LMSSSE <- matrix(,nrow=kfold,ncol=length(targets),byrow =TRUE,dimnames = list(1:kfold,targets))
  # Initialize Rsquared holders
  PLSSSR <- matrix(,nrow=kfold,ncol=length(uniqueTarget),byrow =TRUE,dimnames = list(1:kfold,uniqueTarget))
  PLSSST <- matrix(,nrow=kfold,ncol=length(uniqueTarget),byrow =TRUE,dimnames = list(1:kfold,uniqueTarget))
  PLSRsquared <- matrix(,nrow=reps,ncol=length(uniqueTarget),byrow =TRUE,dimnames = list(1:reps,uniqueTarget))
  # Initialize predMAPE
  PLSSAPE <- matrix(,nrow=kfold,ncol=length(targets),byrow =TRUE,dimnames = list(1:kfold,targets))
  PLSMAPE <- matrix(,nrow=reps,ncol=length(targets),byrow =TRUE,dimnames = list(1:reps,targets))
  LMMAPE <- matrix(,nrow=reps,ncol=length(targets),byrow =TRUE,dimnames = list(1:reps,targets))
  LMSAPE <- matrix(,nrow=kfold,ncol=length(targets),byrow =TRUE,dimnames = list(1:kfold,targets))
  # Initialize predMAD
  PLSSAD <- matrix(,nrow=kfold,ncol=length(targets),byrow =TRUE,dimnames = list(1:kfold,targets))
  PLSMAD <- matrix(,nrow=reps,ncol=length(targets),byrow =TRUE,dimnames = list(1:reps,targets))
  LMMAD <- matrix(,nrow=reps,ncol=length(targets),byrow =TRUE,dimnames = list(1:reps,targets))
  LMSAD <- matrix(,nrow=kfold,ncol=length(targets),byrow =TRUE,dimnames = list(1:kfold,targets))
  
  # Perform repetitions
  for (x in 1:reps) {
    
    #Randomly shuffle the data
    dataset <- dataset[sample(nrow(dataset)),]
    
    # Extract the target and non-target variables for Linear Model
    independentMatrix <- dataset[,sources]
    dependentMatrix <- dataset[,targets]
    
    #Create 10 equally size folds
    folds <- cut(seq(1,nrow(dataset)),breaks=kfold,labels=FALSE)
    
   #Perform 10 fold cross validation
    for(i in 1:kfold){
      #Segment your data by fold using the which() function 
      testIndexes <- which(folds==i,arr.ind=TRUE)
      testingData <- dataset[testIndexes, ]
      trainingData <- dataset[-testIndexes, ]
      indepTestData <- independentMatrix[testIndexes, ]
      indepTrainData <- independentMatrix[-testIndexes, ]
      depTestData <- dependentMatrix[testIndexes, ]
      depTrainData <- dependentMatrix[-testIndexes, ]
      
      #PLS training model
      testHolder <- PLSpredict(trainingData, testingData ,smMatrix, mmMatrix, maxIt, stopCriterion)
      
      #PLS test model for collecting actuals
      PLS_actual_model <- simplePLS(testingData,smMatrix,mmMatrix,maxIt,stopCriterion)
      
      #Initialize PLS residuals and actuals holder matrices
      PLSactuals <- testHolder$testData[,targets]
      PLS_composite_actuals <- as.matrix(PLS_actual_model$fscores[,uniqueTarget])
      PLS_composite_predictions <- as.matrix(testHolder$compositeScores[,uniqueTarget])
      PLSresiduals <- testHolder$residuals[,targets]
      
      #Initialize lm residuals and actuals holder matrices
      lmprediction <- matrix(,nrow=nrow(depTestData),ncol=length(targets),byrow =TRUE,dimnames = list(1:nrow(depTestData),targets))
      lmresidual <- matrix(,nrow=nrow(depTestData),ncol=length(targets),byrow =TRUE,dimnames = list(1:nrow(depTestData),targets))
      lmactual <- matrix(,nrow=nrow(depTestData),ncol=length(targets),byrow =TRUE,dimnames = list(1:nrow(depTestData),targets))
      
      #LM Models
      for(l in 1:length(targets)){
        trainLM <- lm(depTrainData[,l] ~ ., indepTrainData)
        lmprediction[,l] <- predict(trainLM, newdata = indepTestData)
        lmresidual[,l] <- lmprediction[,l] - depTestData[, l]
        lmactual[,l] <- depTestData[, l]
      }
      
      #iterate over no of uniquetargets (endogenous constructs) to calculate oos Rsquared
      for(j in 1:length(uniqueTarget)){
        #Calculate SST and SSR
        PLSSST[i,j] <- sum((PLS_composite_actuals[,j] - mean(PLS_composite_actuals[,j]))^2)
        PLSSSR[i,j] <- sum((PLS_composite_actuals[,j] - PLS_composite_predictions[,j])^2)
        PLSRsquared[x,j] <- 1 - (sum(PLSSSR[,j])/sum(PLSSST[,j]))
      }
      
      #Iterate over no of targets
      for(j in 1:length(targets)){
        
        #Calculate SMSE
        PLSSSE[i,j] <- sum(PLSresiduals[,j]^2) 
        LMSSSE[i,j] <- sum(lmresidual[,j]^2)
        #Calculate SAPE
        PLSSAPE[i,j] <- sum((abs(PLSresiduals[,j]/PLSactuals[,j])))
        #PLSSAPE[i,j] <- sum((abs((mean(testHolder$testData[,j]) - testHolder$predictedMeasurements[,j])/mean(testHolder$testData[,j]))))
        #PLSSAPE[i,j] <- sum((abs((testHolder$residuals[,j])/mean(testHolder$testData[,j]))))
        LMSAPE[i,j] <- sum((abs(lmresidual[,j]/lmactual[,j])))
        #Calculate SAD
        PLSSAD[i,j] <- sum(abs(PLSresiduals[,j] - mean(PLSresiduals[,j])))
        LMSAD[i,j] <- sum(abs(lmresidual[,j] - mean(lmresidual[,j])))
      }
    }
    
    #Final calculations 
    denom <- nrow(dataset)
    for (k in 1:length(targets)) {
      LMRMSE[x,k] <- sqrt((sum(LMSSSE[,k]))/denom)
      PLSRMSE[x,k] <- sqrt((sum(PLSSSE[,k]))/denom)
      LMMAPE[x,k] <- 100*(sum(LMSAPE[,k])/denom)
      PLSMAPE[x,k] <- 100*(sum(PLSSAPE[,k])/denom)
      LMMAD[x,k] <- sum(LMSAD[,k])/denom
      PLSMAD[x,k] <- sum(PLSSAD[,k])/denom
    }
  }
  
  validateResults <- list(PLSRMSE = PLSRMSE, 
                          PLSMAPE = PLSMAPE,
                          PLSMAD = PLSMAD,
                          PLSR2 = PLSRsquared,
                          LMRMSE = LMRMSE,
                          LMMAPE = LMMAPE,
                          LMMAD = LMMAD)
  return(validateResults)
  
}