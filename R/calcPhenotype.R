#ATTENTION: all you need to do to run this script is to determine parameters of the calcPhenotype() function (starts at line 3). 
#__________________________________________________________________________________________________________________________________
#Determine parameters for calcPhenotype() function.
#__________________________________________________________________________________________________________________________________
#Read in training data for GDSC (expression and response)
#_______________________________________________________
#Read in GDSC training expression data. rownames() are genes and colnames() are samples (cell lines/cosmic ids).
trainingExprData=as.matrix(read.table('GDSC_Expression_Matrix.txt', header=TRUE, row.names=1))
dim(trainingExprData) #17419 257 
#Read in GDSC training response data. rownames() are samples (cell lines, cosmic ids), colnames() are drugs.
trainingPtype=as.matrix(read.table('GDSC_Response_Matrix.txt', header=TRUE, row.names=1))
dim(trainingPtype) #257 198
#R puts an X in front of each colname() of expression data :( Fix that here so that cell line data is the same for both matrices.
colnames(trainingExprData)<-rownames(trainingPtype)
#IMPORTANT note: here I do e^IC50 since the IC50s are actual ln values/log transformed already, and the calcPhenotype function Paul has will do a power transformation (I assumed it would be better to not have both transformations)
trainingPtype<-exp(trainingPtype) 

#Or read in training data for CTRP (expression and response)
#_______________________________________________________
#Read in CTRP training expression data. rownames() are genes and colnames() are samples (cell lines/cosmic ids).
trainingExprData=as.matrix(read.table('CTRP_Expression_Matrix.txt', header=TRUE, row.names=1))
dim(trainingExprData) #54356 1076
#Read in CTRP training response data. rownames() are samples (cell lines, cosmic ids), colnames() are drugs.
trainingPtype=as.matrix(read.table('CTRP_Response_Matrix.txt', header=TRUE, row.names=1))
dim(trainingPtype) #1076 545 

#Or read in training data for PRISM (expression and response)
#_______________________________________________________
trainingExprData=readRDS("CCLE_rpkm.rds")
dim(trainingExprData) #56202 1019
trainingPtype=readRDS("prismDrugMat.rds")
dim(t(trainingPtype)) #1419 481 
overlap_cellLines<-intersect(colnames(trainingExprData), rownames(trainingPtype))
trainingExprData=trainingExprData[,overlap_cellLines]
dim(trainingExprData) #56202 473
trainingPtype=trainingPtype[overlap_cellLines,]
dim(trainingPtype) #473 1419

#Read in testing data as a matrix with rownames() as genes and colnames() as samples.
testExprData=as.matrix(read.table('prostate test data.txt', header=TRUE, row.names=1))

#batchCorrect options: "eb" for ComBat, "qn" for quantiles normalization, "standardize", or "none"
#"eb" is good to use when you use microarray training data to build models on microarray testing data.
#"standardize is good to use when you use microarray training data to build models on RNA-seq testing data (this is what Paul used in the 2017 IDWAS paper that used GDSC microarray to impute in TCGA RNA-Seq data, see methods section of that paper for rationale)
batchCorrect<-"eb"

#Determine whether or not to power transform the phenotype data.
#Default is TRUE.
powerTransformPhenotype<-TRUE

#Determine percentage of low varying genes to remove.
#Default is 0.2 (seemingly arbitrary).
removeLowVaryingGenes<-0.2

#Determine method to remove low varying genes.
#Options are 'homogenizeData' and 'rawData'
#homogenizeData is likely better if there is ComBat batch correction, raw data was used in the 2017 IDWAS paper that used GDSC microarray to impute in TCGA RNA-Seq data.
removeLowVaringGenesFrom<-"homogenizeData"

#Determine the minimum number of training samples required to train on.
#Note: this shouldn't be an issue if you train using GDSC or CTRP because there are many samples in both training datasets.
#10, I believe, is arbitrary and testing could be done to get a better number.
minNumSamples=10

#Determine how you would like to deal with duplicate gene IDs.
#Sometimes based on how you clean the data, there shouldn't be any duplicates to deal with.
#Options are -1 for ask user, 1 for summarize by mean, and 2 for disregard duplicates
selection<- 1

#Determine if you'd like to print outputs.
#Default is TRUE.
printOutput=TRUE

#Indicate whether or not you'd like to use PCA for feature/gene reduction. Options are 'TRUE' and 'FALSE'.
#Note: If you indicate 'report_pca=TRUE' you need to also indicate 'pca=TRUE'
pca=FALSE

#Indicate whether you want to output the principal components. Options are 'TRUE' and 'FALSE'.
report_pca=FALSE

#Indicate if you want to convert your training data to TPM, which is recommended if your testing data is measured in TPM. This would be useful, for example, if you wish to use CTRP training data, which is measured in RPKM.
#Options are 'TRUE' and 'FALSE'.
tpm=FALSE

#Indicate if you want correlation coefficients for biomarker discovery. These are the correlations between a given gene of interest across all samples vs. a given drug response across samples.
#These correlations can be ranked to obtain a ranked correlation to determine highly correlated drug-gene associations.
cc=FALSE

#Indicate whether or not you want to output the R^2 values for the data you train on from true and predicted values.
#These values represent the percentage in which the optimal model accounts for the variance in the training data.
#Options are 'TRUE' and 'FALSE'.
rsq=FALSE
#__________________________________________________________________________________________________________________________________
#Set seed.
#__________________________________________________________________________________________________________________________________
set.seed(12345)
#__________________________________________________________________________________________________________________________________
#Load libraries.
#__________________________________________________________________________________________________________________________________
library('sva')
library('preprocessCore')
library('stringr')
library('org.Hs.eg.db')
library('biomaRt')
library('EnsDb.Hsapiens.v75')
library('downloader')
library('car')
library('genefilter')
library('tidyverse')
library('ensembldb')
library('readxl')
library('illuminaHumanv4.db')
library('glmnet')
library('gdata')
library('pls')
#__________________________________________________________________________________________________________________________________
#Read in the calcPhenotype() function as well as the functions that are called in calcPhenotype()
#These functions include doVariableSelection(), homogenizeData(), and summarizeGenesByMean()
#__________________________________________________________________________________________________________________________________
doVariableSelection <- function(exprMat, removeLowVaryingGenes=.2)
{
  vars <- apply(exprMat, 1, var)
  return(order(vars, decreasing=TRUE)[seq(1:as.integer(nrow(exprMat)*(1-removeLowVaryingGenes)))])
}

homogenizeData<-function (testExprMat, trainExprMat, batchCorrect = "eb", selection = -1, printOutput = TRUE)
{
  
  #Check the batchCorrect parameter
  if (!(batchCorrect %in% c("eb", "qn", "none",
                            "rank", "rank_then_eb", "standardize")))
    stop("\"batchCorrect\" must be one of \"eb\", \"qn\", \"rank\", \"rank_then_eb\", \"standardize\" or \"none\"")
  
  #Check if both row and column names have been specified.
  if (is.null(rownames(trainExprMat)) || is.null(rownames(testExprMat))) {
    stop("ERROR: Gene identifiers must be specified as \"rownames()\" on both training and test expression matrices. Both matices must have the same type of gene identifiers.")
  }
  
  #Check that some of the row names overlap between both datasets (print an error if none overlap)
  if (sum(rownames(trainExprMat) %in% rownames(testExprMat)) == 0) {
    stop("ERROR: The rownames() of the supplied expression matrices do not match. Note that these are case-sensitive.")
  }
  else {
    if (printOutput)
      cat(paste("\n", sum(rownames(trainExprMat) %in% rownames(testExprMat)), " gene identifiers overlap between the supplied expression matrices... \n", paste = ""))
  }
  
  #If there are duplicate gene names, give the option of removing them or summarizing them by their mean.
  if ((sum(duplicated(rownames(trainExprMat))) > 0) || sum(sum(duplicated(rownames(testExprMat))) > 0)) {
    if (selection == -1) {
      cat("\nExpression matrix contain duplicated gene identifiers (i.e. duplicate rownames()), how would you like to proceed:")
      cat("\n1. Summarize duplicated gene ids by their mean value (acceptable in most cases)")
      cat("\n2. Disguard all duplicated genes (recommended if unsure)")
      cat("\n3. Abort (if you want to deal with duplicate genes ids manually)\n")
    }
    while (is.na(selection) | selection <= 0 | selection > 3) {
      selection <- readline("Selection: ")
      selection <- ifelse(grepl("[^1-3.]", selection), -1, as.numeric(selection))
    }
    
    cat("\n")
    
    if (selection == 1) #Summarize duplicates by their mean.
    {
      if ((sum(duplicated(rownames(trainExprMat))) > 0)) {
        trainExprMat <- summarizeGenesByMean(trainExprMat)
      }
      if ((sum(duplicated(rownames(testExprMat))) > 0)) {
        testExprMat <- summarizeGenesByMean(testExprMat)
      }
    }
    else if (selection == 2) #Disguard all duplicated genes.
    {
      if ((sum(duplicated(rownames(trainExprMat))) > 0)) {
        keepGenes <- names(which(table(rownames(trainExprMat)) == 1))
        trainExprMat <- trainExprMat[keepGenes, ]
      }
      if ((sum(duplicated(rownames(testExprMat))) > 0)) {
        keepGenes <- names(which(table(rownames(testExprMat)) == 1))
        testExprMat <- testExprMat[keepGenes, ]
      }
    }
    else {
      stop("Exectution Aborted!")
    }
  }
  
  #Subset and order gene ids on the expression matrices.
  commonGenesIds <- rownames(trainExprMat)[rownames(trainExprMat) %in%
                                             rownames(testExprMat)]
  trainExprMat <- trainExprMat[commonGenesIds, ]
  testExprMat <- testExprMat[commonGenesIds, ]
  
  #Subset and order the 2 expression matrices.
  if (batchCorrect == "eb") {
    #Subset to common genes and batch correct using ComBat.
    dataMat <- cbind(trainExprMat, testExprMat)
    mod <- data.frame(`(Intercept)` = rep(1, ncol(dataMat)))
    rownames(mod) <- colnames(dataMat)
    whichbatch <- as.factor(c(rep("train", ncol(trainExprMat)),
                              rep("test", ncol(testExprMat))))
    
    # Filter out genes with low variances to make sure comBat run correctly
    dataMat <- cbind(trainExprMat, testExprMat)
    gene_vars = apply(dataMat, 1, var)
    dataMat = dataMat[-which(gene_vars <= 1e-3),]
    
    combatout <- ComBat(dataMat, whichbatch, mod = mod)
    return(list(train = combatout[, whichbatch == "train"],
                test = combatout[, whichbatch == "test"], selection = selection))
  }
  else if (batchCorrect == "standardize") #Standardize to mean 0 and variance 1 in each dataset using a non EB based approach.
  {
    for (i in 1:nrow(trainExprMat)) {
      row <- trainExprMat[i, ]
      trainExprMat[i, ] <- ((row - mean(row))/sd(row))
    }
    for (i in 1:nrow(testExprMat)) {
      row <- testExprMat[i, ]
      testExprMat[i, ] <- ((row - mean(row))/sd(row))
    }
    return(list(train = trainExprMat, test = testExprMat,
                selection = selection))
  }
  else if (batchCorrect == "rank") #The random-rank transform approach, that may be better when applying models to RNA-seq data.
  {
    for (i in 1:nrow(trainExprMat)) {
      trainExprMat[i, ] <- rank(trainExprMat[i, ], ties.method = "random")
    }
    for (i in 1:nrow(testExprMat)) {
      testExprMat[i, ] <- rank(testExprMat[i, ], ties.method = "random")
    }
    return(list(train = trainExprMat, test = testExprMat,
                selection = selection))
  }
  else if (batchCorrect == "rank_then_eb") #Rank-transform the RNAseq data, then apply ComBat
  {
    #First, rank transform the RNA-seq data.
    for (i in 1:nrow(testExprMat)) {
      testExprMat[i, ] <- rank(testExprMat[i, ], ties.method = "random")
    }
    #Subset to common genes and batch correct using ComBat.
    dataMat <- cbind(trainExprMat, testExprMat)
    mod <- data.frame(`(Intercept)` = rep(1, ncol(dataMat)))
    rownames(mod) <- colnames(dataMat)
    whichbatch <- as.factor(c(rep("train", ncol(trainExprMat)),
                              rep("test", ncol(testExprMat))))
    combatout <- ComBat(dataMat, whichbatch, mod = mod)
    return(list(train = combatout[, whichbatch == "train"],
                test = combatout[, whichbatch == "test"], selection = selection))
  }
  else if (batchCorrect == "qn")
  {
    dataMat <- cbind(trainExprMat, testExprMat)
    dataMatNorm <- normalize.quantiles(dataMat)
    whichbatch <- as.factor(c(rep("train", ncol(trainExprMat)),
                              rep("test", ncol(testExprMat))))
    return(list(train = dataMatNorm[, whichbatch == "train"],
                test = dataMatNorm[, whichbatch == "test"],
                selection = selection))
  }
  else {
    return(list(train = trainExprMat, test = testExprMat,
                selection = selection))
  }
}
summarizeGenesByMean <- function(exprMat)
{
  geneIds <- rownames(exprMat)
  t <- table(geneIds) #How many times each gene name is duplicated.
  allNumDups <- unique(t)
  allNumDups <- allNumDups[-which(allNumDups == 1)]
  
  #Create a *new* gene expression matrix with everything in the correct order....
  #Start by just adding stuff that isn't duplicated
  exprMatUnique <- exprMat[which(geneIds %in% names(t[t == 1])), ]
  gnamesUnique <- geneIds[which(geneIds %in% names(t[t == 1]))]
  
  #Add all the duplicated genes to the bottom of "exprMatUniqueHuman", summarizing as you go
  for(numDups in allNumDups)
  {
    geneList <- names(which(t == numDups))
    
    for(i in 1:length(geneList))
    {
      exprMatUnique <- rbind(exprMatUnique, colMeans(exprMat[which(geneIds == geneList[i]), ]))
      gnamesUnique <- c(gnamesUnique, geneList[i])
      # print(i)
    }
  }
  
  if(class(exprMatUnique) == "numeric")
  {
    exprMatUnique <- matrix(exprMatUnique, ncol=1)
  }
  
  rownames(exprMatUnique) <- gnamesUnique
  return(exprMatUnique)
}
calcPhenotype<-function (trainingExprData,
                         trainingPtype,
                         testExprData,
                         batchCorrect,
                         powerTransformPhenotype=TRUE,
                         removeLowVaryingGenes=0.2,
                         minNumSamples,
                         selection=1,
                         printOutput,
                         removeLowVaringGenesFrom,
                         pca=FALSE, 
                         report_pca=FALSE,
                         rsq=FALSE,
                         cc=FALSE,
                         tpm=FALSE)
{ 
  
  #Initiate empty lists for each data type you'd like to collect.
  #_______________________________________________________________
  DrugPredictions<-list() #Collects drug predictions.
  rsqs<-list() #Collects R^2 values.
  cors<-list() #Collects correlation coefficient for each gene across all samples vs. each drug across all samples.
  
  drugs<-colnames(trainingPtype) #Store all the possible drugs in a vector.
  #drugs=drugs[387]
  
  #Check the supplied data and parameters.
  #_______________________________________________________________
  if (class(testExprData) != "matrix")
    stop("\nERROR: \"testExprData\" must be a matrix.")
  if (class(trainingExprData) != "matrix")
    stop("\nERROR: \"trainingExprData\" must be a matrix.")
  if (class(trainingPtype) != "matrix")
    stop("\nERROR: \"trainingPtype\" must be a matrix.")
  
  if (report_pca)
    if (pca == FALSE)
      stop("\nERROR: pca must be TRUE if report_pca is TRUE")
  
  if (pca)
    if (cc)
      stop("\nERROR: pca must be FALSE if cc is TRUE")
  
  #Make sure training samples are equivalent in both matrices.
  if (colnames(trainingExprData) != rownames(trainingPtype))
    stop("\nERROR: Samples in training matrices must be of equivalent")
  
  #Check if an adequate number of training and test samples have been supplied.
  #_______________________________________________________________
  if ((nrow(trainingExprData) < minNumSamples) || (nrow(testExprData) < minNumSamples)) {
    stop(paste("\nThere are less than", minNumSamples, "samples in your test or training set. It is strongly recommended that you use larger numbers of samples in order to (a) correct for batch effects and (b) fit a reliable model. To supress this message, change the \"minNumSamples\" parameter to this function."))
  }
  
  #If the testing data is in TPM and the training data is in RPKM, convert training data (which is in RPKM) to TPM.
  #_______________________________________________________________
  if (tpm){
    tpm_data <- matrix(0, nrow(trainingExprData), ncol=ncol(trainingExprData), dimnames=trainingExprData)
    for (a in 1:ncol(tpm_data)){ #a represents each gene
      for (b in 1:nrow(tpm_data)){ #b represents each sample
        cell <-trainingExprData[b,a] #Patient b's expression of gene a.
        col_sum <-sum(trainingExprData[,a]) #Sum of that gene's expression across all patients.
        tpm_data[b,a] <- (cell/col_sum) * 10^6
      }
    }
    trainingExprData <- tpm_data
  }
  
  #Get the homogenized data.
  #_______________________________________________________________
  homData <- homogenizeData(testExprData, trainingExprData, batchCorrect, selection, printOutput)
  
  #Remove low varying genes.
  #_______________________________________________________________
  #Do variable selection if specified. By default, we remove 20% of least varying genes from the homogenized dataset.
  #We can also remove the intersection of the lowest 20% from both training and test sets (for the gene ids remaining in the homogenized data).
  #Otherwise, keep all genes.
  
  #Check batchCorrect parameter.
  if (!(removeLowVaringGenesFrom %in% c("homogenizeData", "rawData"))) {
    stop("\nremoveLowVaringGenesFrom\" must be one of \"homogenizeData\", \"rawData\"")
  }
  
  keepRows <- seq(1:nrow(homData$train)) #By default we will keep all the genes.
  if (removeLowVaryingGenes > 0 && removeLowVaryingGenes < 1) { #If the proportion of variability to keep is between 0 and 1.
    if (removeLowVaringGenesFrom == "homogenizeData") { #If you're filtering based on homogenized data.
      keepRows <- doVariableSelection(cbind(homData$test, homData$train), removeLowVaryingGenes = removeLowVaryingGenes)
      
      numberGenesRemoved <- nrow(homData$test) - length(keepRows)
      if (printOutput) cat(paste("\n", numberGenesRemoved, "low variabilty genes filtered."));
    }
    else if (removeLowVaringGenesFrom == "rawData") { #If we are filtering based on the raw data i.e. the intersection of the things filtered from both datasets.
      evaluabeGenes <- rownames(homData$test)
      keepRowsTrain <- doVariableSelection(trainingExprData[evaluabeGenes,], removeLowVaryingGenes = removeLowVaryingGenes)
      keepRowsTest <- doVariableSelection(testExprData[evaluabeGenes,], removeLowVaryingGenes = removeLowVaryingGenes)
      keepRows <- intersect(keepRowsTrain, keepRowsTest)
      numberGenesRemoved <- nrow(homData$test) - length(keepRows)
      if (printOutput)
        cat(paste("\n", numberGenesRemoved, "low variabilty genes filtered."));
    }
  }
  
  #Predict for each drug.
  #_______________________________________________________________
  for(a in 1:length(drugs)){ #For each drug...
    
    #Set up the trainingPtype and trainingExprData. Only use cell lines for which you have response data for.
    #_______________________________________________________________
    trainingPtype2<-trainingPtype[,a] #Obtain the response data for the compound of interest. Must be a numeric vector.
    NonNAindex <- which(!is.na(trainingPtype2)) #Get the indices of the non NAs. You only want the cell lines/cosmic ids that you have drug info for.
    
    samps<-rownames(trainingPtype)[NonNAindex]
    #Index for the drug of interest.
    #Cell lines you will use for training because you have expression and response data for it,
    
    if (length(samps) == 1){
      drugs = drugs[-a]
      cat(paste("\n", drugs[a], "is skipped due to insufficient cell lines to fit the model."))
      next
    }else{
      
      trainingPtype4<-as.numeric(trainingPtype2[NonNAindex])
      
      #PowerTransform the phenotype if specified.
      #_______________________________________________________________
      offset = 0
      if (powerTransformPhenotype){
        if (min(trainingPtype4) < 0){ #All numbers must be positive for a powertransform to work, so make them positive.
          offset <- -min(trainingPtype4) + 1
          trainingPtype4 <- trainingPtype4 + offset
        }
        transForm <- powerTransform(trainingPtype4)[[6]]
        trainingPtype4 <- trainingPtype4^transForm
      }
      
      #Create the ridge regression model on the training data using pca (keeping 2 components).
      #_______________________________________________________________
      if (pca){
        
        #There are many ways for pca in R...here I will use PCR (principal component regression).
        
        #Now use pcr to predict for testing data.
        #_______________________________________________________________
        if (printOutput) cat("\nCalculating predicted phenotype using pca...")
        
        train_x<-(t(homData$train)[samps,keepRows]) #samps represent the cell lines that have been filtered, keepRows represents the genes.
        
        #train_x<-t(homData$train[keepRows,NonNAindex])
        train_y<-trainingPtype4
        data<-as.data.frame(cbind(train_x, train_y))
        
        #Check to make sure you have enough training samples for the drug's model.
        #_______________________________________________________________
        pcrmodel<-try(pcr(train_y~., data=data, validation='CV', ncomp=2))
        if (class(pcrmodel) == 'try-error'){
          drugs = drugs[-a]
          cat(paste("\n", drugs[a], "is skipped due to insufficient cell lines to fit the model."))
          next
        } else{
          #if(printOutput) cat("\nCalculating predicted phenotype...")
          preds<-predict(pcr_model, t(homData$test)[,keepRows], ncomp=2)
          
          #You can compute an R^2 value for the data you train on from true and predicted values.
          #The rsq value represents the percentage in which the optimal model accounts for the variance in the training data.
          #_______________________________________________________________
          if (rsq){
            
            data<-t(homData$train[keepRows,NonNAindex])
            
            if (dim(data)[1] < 4){ #The code will result in an error if you have 3 samples (which is enough for the model fitting but not when you do a 70/30% split)...
              cat(paste("\n", drugs[a], 'is skipped for R^2 analysis'))
            }else{
              data<-(cbind(data, trainingPtype4))
              dt<-sort(sample(nrow(data), nrow(data)*.7)) #sample() randomly picks 70% of rows/samples from the dataset. It samples without replacement.
              
              #Prepare the training data (70% of original training data)
              train_x<-data[dt,]
              ncol<-dim(train_x)[2]
              train_y<-train_x[,ncol]
              train_x<-train_x[,-ncol]
              #Prepare the testing data (30% of original training data)
              test_x<-data[-dt,]
              ncol<-dim(test_x)[2]
              test_y<-test_x[,ncol]
              test_x<-test_x[,-ncol]
              
              data<-as.data.frame(cbind(train_x, train_y))
              pcr_model<-pcr(train_y~., data=data, validation='CV', ncomp=2) #Set validation argument to CV.
              pcr_pred<-predict(pcr_model, test_x, ncomp=2)
              
              if (printOutput) cat("\nCalculating R^2...")
              sst<-sum((pcr_pred - mean(test_y))^2) #Compute the sum of squares total.
              sse<-sum((pcr_pred - test_y))^2 #Compute the sum of squares error.
              rsq_value<-1 - sse/sst #Compute the rsq value.
            }
          }
        }
        
        if (report_pca){
          if (printOutput) cat("\nObtaining principal components...")
          pcs<-coef(pcr_model, comps = c(1,2)) #comps: numeric, which components to return.
          dir.create("./calcPhenotype_Output")
          path<-paste('./calcPhenotype_Output/', drugs[a], '.RData', sep="")
          save(pcs, file=path)
        }
        
      } else {
        
        #Create the ridge regression model on our training data to predict for our actual testing data without pca.
        #_______________________________________________________________
        if(printOutput) cat("\nFitting Ridge Regression model...");
        
        expression<-(t(homData$train)[samps,keepRows]) #samps represent the cell lines that have been filtered, keepRows represents the genes.
        
        rrModel<-glmnet(expression, trainingPtype4, family='gaussian', alpha=0)
        
        #Check to make sure you have enough training samples for the drug's model.
        #_______________________________________________________________
        cv_fit<-try(cv.glmnet(expression, trainingPtype4, alpha=0), silent = TRUE)
        if (class(cv_fit) == 'try-error'){
          drugs = drugs[-a]
          cat(paste("\n", drugs[a], "is skipped due to insufficient cell lines to fit the model."))
          next
        } else{
          if(printOutput) cat("\nCalculating predicted phenotype...")
          preds<-predict(rrModel, newx=t(homData$test)[,keepRows], s=cv_fit$lambda.min)
          
          #You can compute an R^2 value for the data you train on from true and predicted values.
          #The R^2 value represents the percentage in which the optimal model accounts for the variance in the training data.
          #_______________________________________________________________
          if(rsq){
            
            data<-t(homData$train[keepRows,NonNAindex])
            
            if (dim(data)[1] < 4){ #The code will result in an error if you have 3 samples...
              cat(paste("\n", drugs[a], 'is skipped for R^2 analysis'))
            }else{
              data<-(cbind(data, trainingPtype4))
              dt<-sort(sample(nrow(data), nrow(data)*.7)) #sample() randomly picks 70% of rows/samples from the dataset. It samples without replacement.
              
              #Prepare the training data (70% of original training data)
              train_x<-data[dt,]
              ncol<-dim(train_x)[2]
              train_y<-train_x[,ncol]
              train_x<-train_x[,-ncol]
              #Prepare the testing data (30% of original training data)
              test_x<-data[-dt,]
              ncol<-dim(test_x)[2]
              test_y<-test_x[,ncol]
              test_x<-test_x[,-ncol]
              
              rrModel<-glmnet(train_x, train_y, family='gaussian', alpha=0)
              cv_fit<-cv.glmnet(train_x, train_y, alpha=0)
              pred<-predict(rrModel, newx=test_x, s=cv_fit$lambda.min)
              
              if(printOutput) cat("\nCalculating R^2...")
              sst<-sum((pred - mean(test_y))^2) #Compute the sum of squares total.
              sse<-sum((pred - test_y))^2 #Compute the sum of squares error.
              rsq_value<-1 - sse/sst #Compute the rsq value.
            }
          }
        }
      }
      
      #If the response variable was transformed (aka powerTransformPhenotype=TRUE), untransform it here.
      #_______________________________________________________________
      if(powerTransformPhenotype) {
        preds <- preds^(1/transForm)
        preds <- preds - offset
      }
      
      #Find correlation between imputed response and expression of a gene.
      #Each gene is corrected differently; therefore, it may not be ideal to determine weight or percentage or weight in which each gene contributes to the prediction.
      #_______________________________________________________________
      if(cc){ #You can only collect correlations if pca=FALSE!
        if(pca){
          stop('ERROR: pca must equal FALSE in order to compute correlations') #It doesn't make sense to compute correlations when the features have changed from genes to pcs.
        }
        
        if(printOutput) cat("\nCalculating correlation coefficients...") #This is only relevant if you aren't using pca.
        
        cors_vec<-c()
        matrix<-homData$test[keepRows,] #Matrix of genes x cell lines/cosmic ids.
        for(d in 1:nrow(matrix)){ #For each gene...
          cors_vec[d]<-cor(as.vector(matrix[d,]), as.vector(preds)) #Compute correlation coefficient for expression of a given gene across patients vs.imputed values for a given drug across patients.
        }
        #indices<-order(cor, decreasing=TRUE)
        #ordered_cor<-cor[indices] #Order the correlation coefficients from big to small.
        #genes<-rownames(matrix)
        #ordered_genes<-genes[indices]
        #names(ordered_cor)<-ordered_genes
      }
      
      if(printOutput) cat(paste("\nDone making prediction for drug", a, "of", ncol(trainingPtype)))
      
      #Store the data in your lists.
      #_______________________________________________________________
      DrugPredictions[[a]]<-preds
      
      if(rsq){
        rsqs[[a]]<-rsq_value
      }
      
      if(cc){
        cors[[a]]<-cors_vec
      }
    }
  }
  
  #Time to save the data!
  #_______________________________________________________________
  #Save drug prediction data to your home directory as a .txt file.
  names(DrugPredictions)<-drugs
  DrugPredictions_mat<-do.call(cbind, DrugPredictions)
  colnames(DrugPredictions_mat)<-drugs
  rownames(DrugPredictions_mat)<-colnames(testExprData)
  dir.create("./calcPhenotype_Output")
  write.table(DrugPredictions_mat, file="./calcPhenotype_Output/DrugPredictions.txt")
  
  #If rsq=TRUE, save R^2 data.
  if(rsq){
    names(rsqs)<-drugs
    rsqs_mat<-do.call(cbind, rsqs)
    dir.create("./calcPhenotype_Output")
    write.table(rsqs_mat, file="./calcPhenotype_Output/R^2.txt")
  }
  
  #If CC=TRUE, save correlation coefficient data.
  if(cc){
    names(cors)<-drugs
    cor_mat<-do.call(cbind, cors)
    rownames(cor_mat)<-rownames(homData$train[keepRows,NonNAindex])
    colnames(cor_mat)<-drugs
    dir.create("./calcPhenotype_Output")
    write.table(cor_mat, file="./calcPhenotype_Output/cors.txt")
  }
}
#__________________________________________________________________________________________________________________________________
#Run the calcPhenotype() function using the parameters you specified.
#__________________________________________________________________________________________________________________________________
calcPhenotype(trainingExprData=trainingExprData,
              trainingPtype=trainingPtype,
              testExprData=testExprData,
              batchCorrect=batchCorrect,
              powerTransformPhenotype=powerTransformPhenotype,
              removeLowVaryingGenes=removeLowVaryingGenes,
              minNumSamples=minNumSamples,
              selection=selection,
              printOutput=printOutput,
              removeLowVaringGenesFrom=removeLowVaringGenesFrom,
              pca=pca,
              report_pca=report_pca,
              rsq=rsq,
              cc=cc,
              tpm=tpm)
