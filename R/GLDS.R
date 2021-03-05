#'This function performs an iterative matrix completion algorithm to predict drug response for pre-clinical data when there are missing ('NA') values.
#'@param senMat A matrix of drug sensitivity data with missing ('NA') values. rownames() are samples (e.g. cell lines), and colnames() are drugs.
#'@param nPerms The number of iterations that the EM-algorithm (expectation maximization approach)  run. The default is 50, as previous findings recommend 50 iterations (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1050-9)
#'@return A matrix of drug sensitivity scores without missing values. rownames() are samples, and colnames are drugs.
#'@keywords Drug response prediction.
#'@import glmnet
#'@export
completeMatrix <- function(senMat, nPerms=50)
{
  cat("\nNumber of iterations:")
  #To initialize the algorithm, all missing values are first imputed to the median.
  numCellLinesNotScreened <- apply(senMat, 2, function(r)return(sum(is.na(r))))
  numDrugsNotScreened <- apply(senMat, 1, function(r)return(sum(is.na(r))))
  indexNotScreened <- apply(senMat, 2, function(r)return(which(is.na(r))))
  drugMed <- apply(senMat, 2, function(col)return(median(na.omit(as.numeric(col)))))
  hundIc50sImpute <- data.matrix(senMat)
  for(i in 1:ncol(senMat))
  {
    datImpute <- as.numeric(senMat[,i])
    datImpute[is.na(datImpute)] <- drugMed[i]
    hundIc50sImpute[,i] <- datImpute
  }

  #Sort the matrix by least number of missing values to most.
  #The rows (e.g. cell lines) of the senMat matrix must be ordered, ascending by the number of missing values (so the drugs toward the bottom have more missing values
  #then the drugs at the top).
  hundIc50sImputeSort <- hundIc50sImpute[, order(numCellLinesNotScreened[colnames(hundIc50sImpute)])]

  #Calculate/estimate the PCs of this matrix X.
  matrix<-matrix(as.numeric(unlist(hundIc50sImputeSort)),nrow=nrow(hundIc50sImputeSort))
  colnames(matrix)<-colnames(hundIc50sImputeSort)
  rownames(matrix)<-rownames(hundIc50sImputeSort)
  hundIc50sImputeSort<-matrix
  pcs <- prcomp(hundIc50sImputeSort)$x

  #Fit a lasso model for all other drugs, for all samples other than those NOT screened by the drug we are predicting for. By default, we will iterate 50 times.
  #For each missing value of X, fit a lasso regression model using the PCs of all other rows of X as predictors of all other values (both measured and imputed values)
  #for that cell line.
  imputeSortList <- list()
  imputeSortList[[1]] <- hundIc50sImputeSort
  medianDistance <- numeric()
  for(j in 1:nPerms)
    #Repeat this procedure iteratively until the total change in A (the sum of the total difference between each of the elements in X and X') converges.
    #This takes approximately 50 iterations, which is why 50 is the default.
  {
    pb <- txtProgressBar(min = 0, max = ncol(hundIc50sImputeSort), style = 3) #Create progress bar.
    for(i in 1:ncol(hundIc50sImputeSort))
    {
      setTxtProgressBar(pb, i) #Update progress bar.
      pcs <- prcomp(hundIc50sImputeSort)$x # calcualte the PCs of the current matrix #hundIc50sImputeSort.
      indexNotScreenedThisDrug <- indexNotScreened[[colnames(hundIc50sImputeSort)[i]]] #Index of "imputed" cell lines for this drug.

      if(length(indexNotScreenedThisDrug) > 0)
      {
        #Make the training (design matrix) for the non-imputed cell lines for this drug.
        pcMat <- pcs[-indexNotScreenedThisDrug, ]
        dmRaw <-  model.matrix(~pcMat)

        #Calculate lambda using CV.
        #The tuning parameter lambda for the lasso regression is selected using cross-validation (cv.glment)
        FitRaw <- cv.glmnet(dmRaw, hundIc50sImputeSort[-indexNotScreenedThisDrug, i], nfolds=20)
        getLambdasRaw <- FitRaw$lambda.min #The optimal lambda value.

        #Fit the model and extract the co-efficients.
        theModRaw <- glmnet(dmRaw, hundIc50sImputeSort[-indexNotScreenedThisDrug, i], lambda=getLambdasRaw)
        coef(theModRaw)[,1][coef(theModRaw)[,1] != 0]
        betas <- coef(theModRaw)[,1][coef(theModRaw)[,1] != 0][-1]
        intercept <- coef(theModRaw)[,1][coef(theModRaw)[,1] != 0][1]
        names(betas) <- substring(names(betas), 6)

        #Use the model to update the estimates for the imputed values for this drug.
        #Apply the model to yield an updated estimate for each missing value.
        #Repeat this procedure for all missing values of X to yield an updated matrix X'
        if(length(indexNotScreenedThisDrug) == 1)
        {
          predictions <- intercept + apply((t(pcs[indexNotScreenedThisDrug, names(betas)]) * betas), 1, sum); # P redict new values
        }
        else
        {
          predictions <- intercept + apply(t(t(pcs[indexNotScreenedThisDrug, names(betas)]) * betas), 1, sum);
        }

        hundIc50sImputeSort[indexNotScreenedThisDrug, i] <- predictions  # Here is the updated matrix
      }

    }
    close(pb)
    imputeSortList[[j + 1]] <- hundIc50sImputeSort
    medianDistance[j] <- mean(as.numeric(imputeSortList[[j]]) - as.numeric(imputeSortList[[j + 1]])) #Estimate A (the sum of the total difference between each of the elements of X and X').
    cat(paste("\nIteration: ", j, "\n"), "")
    #plot(medianDistance, main=j, ylab="Median Distance from previous imputed matrix")
  }
  cat("\nDone\n")

  return(hundIc50sImputeSort[,colnames(senMat)])
}
#'This function determines drug-gene associations for pre-clinical data.
#'@param drugMat A matrix of drug sensitivity data. rownames() are pre-clinical samples, and colnames() are drug names.
#'@param drugRelationshipList A list, which for each drug, contains a vector of 1's and 0's, indicating whether the drugs are related (e.g. if both drugs were an inhibitor of ERBB2, that position would contain a 1).
#'The order of the drugs should be the same as the order of the drug sensitivity matrix drugMat.
#'@param markerMat A matrix containing the data for which you are looking for an association with drug sensitivity (e.g. a matrix of somatic mutation data). rownames() are marker names (e.g. gene names), and colnames() are samples.
#'@param numCorDrugsExclude The number of highly correlated drugs to remove. The default is 100.
#'#When calculating GLDS, it is recommended to exclude 25% of drugs in your data set, which have a high correlation with the drug of interest, from the calculation of GLDS to prevent removing a drug specific signal.
#'@param minMuts The minimum number of non-zero entries required so that a p-value can be calculated (e.g. how many somatic mutations must be present). The default is 5.
#'@param additionalCovariateMatrix A matrix containing covariates to be fit in the drug biomarker association models. This could be, for example, tissue of origin or cancer type. Columns are sample names. The default is NULL.
#'@export
gldsCorrectedAssoc <- function(drugMat, drugRelationshipList, markerMat, numCorDrugsExclude=100, minMuts=5, additionalCovariateMatrix=NULL)
{
  results_gldsPs <- list()
  results_gldsBetas <- list()
  results_naivePs <- list()
  results_naiveBetas <- list()
  numDrugs <- ncol(drugMat)
  pairCor <- cor(drugMat, method="spearman")
  comNames <- colnames(markerMat)[colnames(markerMat) %in% rownames(drugMat)] #Cell lines for which we have both mutation and drug data....

  if(!is.null(additionalCovariateMatrix)) #If additional co variate matrix was provided, then also subset to those samples.
  {
    comNames <- comNames[comNames %in% rownames(additionalCovariateMatrix)]
  }

  pb <- txtProgressBar(min = 0, max = ncol(drugMat), style = 3) # Create progress bar
  for(i in 1:ncol(drugMat))
  {

    #Estimate the GLDS as the first 10 PCs of this set of negative controls/non-related sets of drugs.
    #These PCs are included as covariates in the linear model used to perform the subsequent corrected cancer gene mutation to drug IC50 association analysis.
    negControlDrugs <- colnames(drugMat)[!as.logical(drugRelationshipList[[i]])]
    pairwiseCorNear <- names(rank(abs(pairCor[, colnames(drugMat)[i]]))[(numDrugs-numCorDrugsExclude):numDrugs]) #Also remove very correlated drugs. Number defined by "numCorDrugsExclude".
    negControlDrugs <- setdiff(negControlDrugs, pairwiseCorNear) #Remove very highly correlated drugs from "negative controls.
    controlPCsAll <- prcomp(drugMat[, negControlDrugs])$x
    controlPCsAllCom <- controlPCsAll[comNames, ]

    #Calculate the P-values and beta values for each marker for this drug, controlling for GLDS and not controlling for GLDS.
    results_gldsPs[[i]] <- numeric()
    results_gldsPs[[i]] <- rep(NA, nrow(markerMat))
    results_gldsBetas[[i]] <- numeric()
    results_gldsBetas[[i]] <- rep(NA, nrow(markerMat))
    results_naivePs[[i]] <- numeric()
    results_naivePs[[i]] <- rep(NA, nrow(markerMat))
    results_naiveBetas[[i]] <- numeric()
    results_naiveBetas[[i]] <- rep(NA, nrow(markerMat))
    for(j in 1:nrow(markerMat)){
      if(sum(markerMat[j, comNames]) > minMuts){ #comNames are the common cell line names.
        if(is.null(additionalCovariateMatrix)) #If no additional covariate have been provided....
        {
          theCoefs <- coef(summary(lm(drugMat[comNames, i]~markerMat[j, comNames])))
          #Fit the linear models using the standard lm() function.
          #Calculate the p-values using the summary() function, which calculates the significance of the model coefficients using t-tests.
          results_naivePs[[i]][j] <- theCoefs[2, 4]
          results_naiveBetas[[i]][j] <- theCoefs[2, 1]

          theCoefs <- coef(summary(lm(drugMat[comNames, i]~markerMat[j, comNames]+controlPCsAllCom[,1:10])))
          results_gldsPs[[i]][j] <- theCoefs[2, 4]
          results_gldsBetas[[i]][j] <- theCoefs[2, 1]
        }
        else #If there are other covariates, include them in the models.
        {
          theCoefs <- coef(summary(lm(drugMat[comNames, i]~markerMat[j, comNames]+additionalCovariateMatrix[comNames,])))
          #Fit the linear models using the standard lm() function.
          #Calculate the p-values using the summary() function, which calculates the significance of the model coefficients using t-tests.
          results_naivePs[[i]][j] <- theCoefs[2, 4]
          results_naiveBetas[[i]][j] <- theCoefs[2, 1]

          theCoefs <- coef(summary(lm(drugMat[comNames, i]~markerMat[j, comNames]+controlPCsAllCom[,1:10]+additionalCovariateMatrix[comNames,])))
          results_gldsPs[[i]][j] <- theCoefs[2, 4]
          results_gldsBetas[[i]][j] <- theCoefs[2, 1]
        }
      }
    }
    #cat(paste(i, " ", sep=""))
    names(results_gldsPs[[i]]) <- rownames(markerMat)
    names(results_naivePs[[i]]) <- rownames(markerMat)
    names(results_gldsBetas[[i]]) <- rownames(markerMat)
    names(results_naiveBetas[[i]]) <- rownames(markerMat)

    setTxtProgressBar(pb, i)#Update progress bar

  }
  close(pb)
  names(results_gldsPs) <- colnames(drugMat)
  names(results_naivePs) <- colnames(drugMat)
  names(results_gldsBetas) <- colnames(drugMat)
  names(results_naiveBetas) <- colnames(drugMat)

  outList <- list(pGlds=results_gldsPs, betaGlds=results_gldsBetas, pNaive=results_naivePs, betaNaive=results_naiveBetas)
  write.table(outList, file='./GLDS_output.txt')
  return(outList)
}