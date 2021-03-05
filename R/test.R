#'This function will test every drug against every CNV or somatic mutation for your cancer type.
#'@param drug_prediction The drug prediction data. Must be a data frame. rownames() are samples, and colnames() are drugs. Make sure that the rownames() are of the same form as the sample names in your cnv or mutation data.
#'e.g. if the rownames() are TCGA barcodes of the form TCGA-##-####-###, make sure your cnv/mutation data also uses samples in the form TCGA-##-####-###
#'@param data The cnv or mutation data. Must be a data frame. If you wish to use cnv data, use the output from map_cnv().
#'If you wish to use mutation data, use the method for downloading mutation data outlined in the vignette; if not, ensure your data file
#'includes the following columns: 'Variant_Classification', 'Hugo_Symbol', 'Tumor_Sample_Barcode'.
#'@param n The number of samples you want CNVs or mutations to be amplified in. The default is 10.
#'@param cnv TRUE or FALSE. Indicate whether or not you would like to test cnv data. If TRUE, you will test cnv data. If FALSE, you will test mutation data.
#'@keywords Test CNV or mutation data to genes.
#'@import org.Hs.eg.db
#'@import gdata
#'@import parallel
#'@import TxDb.Hsapiens.UCSC.hg19.knownGene
#'@import GenomicFeatures
#'@export
test<-function(drug_prediction, data, n=10, cnv){

  #Check parameters.
  #_____________________________________________________________________________
  if (class(drug_prediction) != "data.frame")
    stop("\nERROR: \"drug_prediction\" must be a data frame")
  if (class(data) != "data.frame")
    stop("\nERROR: \"data\" must be a data frame")


  #If TCGA is in my colnames() (as it would if you got cnv data from map_cnv() OR
  #if TCGA is in the column of your mutation data mutation$Tumor_Sample_Barcode, then you have TCGA samples
  #and you want to make sure you only use 01A samples).
  cols<-colnames(data)
  #_____________________________________________________________________________
  if ( (length(grepl('TCGA', cols)) > 0) | (length(grepl('TCGA', data$Tumor_Sample_Barcode)) > 0) ){ #If the patient samples are TCGA barcodes, then obtain the primary tumor/01A samples.
    if (cnv){
      #If cnv=TRUE, then you have cnv data, so proceed with this testing...
      #Find the indices of the 01A/primary tumor samples and the non-duplicates.
      indices <- which(sapply(strsplit(cols, "-"), function(a)a[4]) == "01A")
      #Create a matrix using only the primary tumor patients as columns.
      matrix <- data[, indices] #Make a matrix of rows/genes and columns/TCGA O1A samples.

      #Collect all the patient names/rows of the drug prediction matrix that are 01A primary tumors.
      drugs_01A <- rownames(drug_prediction)[which(sapply(strsplit(rownames(drug_prediction), "-"), function(a)a[4]) == "01A")]
      #Collect the patient IDs of all the patients in drugs_01A (the patients IDS for the 01A samples you have drug prediction data for). These are the 4 digits after the bacode TCGA-## so TCGA-##-????
      drugs_01Aids <- sapply(strsplit(drugs_01A, "-"), function(a)a[3])
      #Index the drug prediction matrix so that its rows/patients are the primary tumor patients.
      matrix2 <- drug_prediction[drugs_01A,]
      #Rename the rows of the drug prediction matrix to the patient id (the 4 digits).

      #Make sure you only use unique ids (sometimes there are duplicates). If there are duplicates, remove.
      indices<-match(unique(drugs_01Aids), drugs_01Aids)
      matrix2<-matrix2[indices,] #Only keep the rows that have unique ids.
      rownames(matrix2) <- drugs_01Aids[indices]

      #Collect all the patient IDs of the patients in columns of matrix
      ids <- sapply(strsplit(colnames(matrix), "-"), function(a)a[3])

      #Collect the patients that are common to both prediction and cnv/mut data (using their 4 digit ids).
      overlapping_ids<-intersect(drugs_01Aids, ids)

      #Index the drug prediction matrix so that its rows only contains patients it has in common with the cnv/mut data.
      indices<-match(overlapping_ids, drugs_01Aids)
      drug_prediction2<-drug_prediction[indices,]

      indices<-match(overlapping_ids, ids)
      matrix3<-as.matrix(matrix[,indices])

      MatCommonPats_amps <- apply(matrix3, 2, function(theCol)return(as.numeric(theCol > 1))) #Apply this function to each column.
      rownames(MatCommonPats_amps) <- rownames(matrix3)

      MatCommonPats_dels <- apply(matrix3, 2, function(theCol)return(as.numeric(theCol < -1)))

      rownames(MatCommonPats_dels) <- rownames(matrix3)

      theFun <- function(j){
        pVals <- numeric()
        betaVal <- numeric()
        if(sum(na.omit(MatCommonPats_amps[j, ])) > n){ #Make sure the gene is amplifed at least 50 times
          for(i in 1:ncol(drug_prediction2)){
            theMod <- coef(summary(lm(drug_prediction2[,i]~MatCommonPats_amps[j, ]))) #Linear regression betwween each amplified CNV (j/row) and each drug (i/column).
            pVals[i] <- theMod[2,4] #P-values for each model.
            betaVal[i] <- theMod[2,1] #Effect size for each model.
          }
          names(pVals) <- colnames(drug_prediction2)
          names(betaVal) <- colnames(drug_prediction2)
        }
        return(list(pVals, betaVal))
      }

      #mclapply is a parallelized vector of lapply.
      #It's structure is mclapply(X, FUN...) where X is a vector and FUN is the function applied to each element of X.
      #Here, the function is applied to each value of X which is j (each CNV).
      allCors <- mclapply(1:nrow(matrix3), theFun) #You may want to change this number based on the number of available cores.
      #'allCors'is a list where each element/drug consists of 2 elements (the p value for that drug's association with the CNV and the beta value).

      #Name each element after each CNV.
      names(allCors) <- rownames(matrix3)

      hasAmps <- apply(MatCommonPats_amps, 1, function(theRow)return(sum(na.omit(theRow)) > n)) #Restrict analysis to CNAs that occur in 50 or more samples.

      allCors_hasAmps <- allCors[hasAmps]

      pVals <- sapply(allCors_hasAmps, function(item)return(item[[1]]))
      betas <- sapply(allCors_hasAmps, function(item)return(item[[2]]))

      #output<-cbind(pVals, betas)

      write.table(pVals, file='./CnvTestOutput_pVals.txt')
      write.table(betas, file='./CnvTestOutput_betas.txt')

      #return((output)) # its going to be difficult to get at causality in a systematic way here....

    }else{
      #Obtain a list for the patients you have data for and initiate empty lists to fill.
      #_______________________________________
      tcgaIds<-data$Tumor_Sample_Barcode
      unique<-unique(tcgaIds)
      mutsListAll <- list() #A list of all the mutation occurring in each sample.
      proteinChangingMutations <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "In_Frame_Del",
                                    "Frame_Shift_Ins", "In_Frame_Ins", "Nonstop_Mutation", "De_novo_Start_OutOfFrame",
                                    "De_novo_Start_InFrame", "Missense", "Read-through", "Indel")
      genesWithProteinChangeList <- list()

      #Fill those lists.
      #_______________________________________
      for(i in 1:length(unique)){
        print(paste(i, "of", length(unique)), sep='')
        indices<-unique[i] == tcgaIds

        #Now for each sample, pull out a list of genes with somatic mutations
        variantType <- data[indices,"Variant_Classification"]
        theGenes <- data[indices,"Hugo_Symbol"]
        names(variantType) <- theGenes
        mutsListAll[[i]] <- variantType #A list of all the mutation occurring in each sample.
        genesWithProteinChangeList[[i]] <- unique(names(variantType[variantType %in% proteinChangingMutations]))
      }

      allMutatedGenes <- unique(names(unlist(mutsListAll))) #All the unique genes mutated across all patients.
      mutationTypes <- table(unlist(mutsListAll)) #Frequency of mutation type across all patients.
      mutsListAll_unlist <- unlist(mutsListAll, recursive=F) #Mutation type and gene association.
      genesWithProteinChangeList_unlist <- unlist(genesWithProteinChangeList, recursive=F) #Genes with protein changing mutations across all patients.

      #From mutsListAll we can then create a matrix indicating if the gene has a coding mutation.
      #_______________________________________
      mutMat <- numeric((length(unique)*length(allMutatedGenes)))
      dim(mutMat) <- c(length(allMutatedGenes), length(unique))
      rownames(mutMat) <- allMutatedGenes
      colnames(mutMat) <- unique

      #Now populate this matrix with the relevant information about what kind of mutation each gene has in each sample.
      #_______________________________________
      for(i in 1:length(unique)){
        print(paste(i, "of", length(unique)), sep='')
        mutMat[genesWithProteinChangeList[[i]], i] <- rep(1, length(genesWithProteinChangeList[[i]]))
      }

      tumorTypeId <- sapply(strsplit(colnames(mutMat), "-", fixed=TRUE), function(l)return(l[4]))

      #Lets remove everything but the "Primary Solid Tumors (i.e. "01")".
      #_______________________________________
      mutMat_only01 <- mutMat[, tumorTypeId == "01A"] #Genes that were mutated in 01A samples.
      theIds <- colnames(mutMat_only01) #Patient ids that were 01A.
      mutIds <- sapply(strsplit(theIds, "-", fixed=T), function(l)return(l[3])) #The TCGA #### part of the barcode.
      colnames(mutMat_only01) <- mutIds

      #Extract the 01a samples from the drug prediction data, i.e. tumor samples.
      #_______________________________________
      all01ASamples <- colnames(drug_prediction)[which(sapply(strsplit(colnames(drug_prediction), ".", fixed=T), function(a)a[4]) == "01A")]
      preds01a <- drug_prediction[, all01ASamples] #Predictions for 01A samples.
      sampIds01a <- sapply(strsplit(all01ASamples, ".", fixed=T), function(l)return(l[3])) #The TCGA #### digit number.
      colnames(preds01a) <- sampIds01a
      inPredAndMutData <- sampIds01a[sampIds01a %in% mutIds] #Samples for which we have both predicted drug response and mutation calls

      #Run the associations between all genes and drugs, for drugs with at least 50 mutations.
      #_______________________________________
      preds01a_filt_ord <- preds01a[, inPredAndMutData] #The preds for the 01A samples we have both prediction and mutation data for.
      mutMat_nodups_ordFilt <- mutMat_only01[, inPredAndMutData]
      commonMuts <- apply(mutMat_nodups_ordFilt, 1, sum)
      if (length(which(commonMuts >= n)) == 0){
        stop((paste("\nERROR: Mutations were not identified in at least", n, "genes. Recommend decreasing the n parameter.", sep=" ")))
      }
      commonlyMutated <- mutMat_nodups_ordFilt[which(commonMuts >= n), ]

      #If there are gene entries with an unknown HUGO ID, remove it.
      #_______________________________________
      if("Unknown" %in% rownames(commonlyMutated)){
        indices<-'Unknown' %in% rownames(commonlyMutated)
        commonlyMutated<-commonlyMutated[-indices,]
      }

      #Get p values and beta values.
      #_______________________________________
      pValList <- list()
      betaValList <- list()
      for(i in 1:nrow(preds01a_filt_ord)){ #For each drug...
        pValList[[i]] <- numeric()
        betaValList[[i]] <- numeric()
        for(j in 1:nrow(commonlyMutated))
        {
          thecoefs <- coef(summary(lm(preds01a_filt_ord[i,]~commonlyMutated[j,])))
          pValList[[i]][[j]] <- thecoefs[2,4]
          betaValList[[i]][[j]] <- thecoefs[2,1]
        }
      }

      #Get the adjusted p-value for each gene-drug combination, pull out the significant associations
      #and create a supplementary table that lists these for "predictable" drugs?.
      #_______________________________________
      sigPs <- list()
      pAdjListCantype <- list()
      for(i in 1:length(pValList))
      {
        names(pValList[[i]]) <- rownames(commonlyMutated)
        names(betaValList[[i]]) <- rownames(commonlyMutated)
        padj <- p.adjust(pValList[[i]], method="BH")
        sigPs[[i]] <- padj[padj < 0.05]
        pAdjListCantype[[i]] <- padj
      }
      names(sigPs) <- rownames(preds01a_filt_ord)
      names(pValList) <- rownames(preds01a_filt_ord)
      names(betaValList) <- rownames(preds01a_filt_ord)
      names(pAdjListCantype) <- rownames(preds01a_filt_ord)

      output<-sort(unlist(pValList))[1:30]

      write.table(output, file='./MutationTestOutput.txt')

      #Print the top associations
      return(output)
    }
    #This code is for when you don't have TCGA barcoded samples (it's similar to above)
  } else {
    if (cnv){
      overlapping_samples<-intersect(rownames(drug_prediction), colnames(data))

      #Index the drug prediction matrix so that it contains the overlapping samples.
      indices<-match(overlapping_samples, rownames(drug_prediction))
      drug_prediction2<-drug_prediction[indices,]
      #Index the cnv.mutation matrix so that it contains the overlapping samples.
      indices<-match(overlapping_samples, colnames(data))
      matrix3<-as.matrix(matrix[,indices])

      MatCommonPats_amps <- apply(matrix3, 2, function(theCol)return(as.numeric(theCol > 1))) #Apply this function to each column.
      rownames(MatCommonPats_amps) <- rownames(matrix3)

      MatCommonPats_dels <- apply(matrix3, 2, function(theCol)return(as.numeric(theCol < -1)))

      rownames(MatCommonPats_dels) <- rownames(matrix3)

      theFun <- function(j){
        pVals <- numeric()
        betaVal <- numeric()

        if(sum(na.omit(MatCommonPats_amps[j, ])) > n){ # Make sure the gene is amplifed at least 50 times
          for(i in 1:ncol(drug_prediction2)){
            theMod <- coef(summary(lm(drug_prediction2[,i]~MatCommonPats_amps[j, ]))) #Linear regression betwween each amplified CNV (j/row) and each drug (i/column).
            pVals[i] <- theMod[2,4] #P-values for each model.
            betaVal[i] <- theMod[2,1] #Effect size for each model.
          }
          names(pVals) <- colnames(drug_prediction2)
          names(betaVal) <- colnames(drug_prediction2)
        }
        return(list(pVals, betaVal))
      }

      #mclapply is a parallelized vector of lapply.
      #It's structure is mclapply(X, FUN...) where X is a vector and FUN is the function applied to each element of X.
      #Here, the function is applied to each value of X which is j (each CNV).
      allCors <- mclapply(1:nrow(matrix3), theFun) #You may want to change this number based on the number of available cores.
      #'allCors'is a list where each element/drug consists of 2 elements (the p value for that drug's association with the CNV and the beta value).

      #Name each element after each CNV.
      names(allCors) <- rownames(matrix3)

      hasAmps <- apply(MatCommonPats_amps, 1, function(theRow)return(sum(na.omit(theRow)) > n)) #Restrict analysis to CNAs that occur in 50 or more samples.

      allCors_hasAmps <- allCors[hasAmps]

      pVals <- sapply(allCors_hasAmps, function(item)return(item[[1]]))
      betas <- sapply(allCors_hasAmps, function(item)return(item[[2]]))

      write.table(pVals, file='./CnvTestOutput.txt')

      return((pVals)) # its going to be difficult to get at causality in a systematic way here....

    } else {

      #Obtain a list for the patients you have data for and initiate empty lists to fill.
      #_______________________________________
      sampIds<-data$Tumor_Sample_Barcode
      unique<-unique(sampIds)
      mutsListAll <- list() #A list of all the mutation occurring in each sample.
      proteinChangingMutations <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "In_Frame_Del",
                                    "Frame_Shift_Ins", "In_Frame_Ins", "Nonstop_Mutation", "De_novo_Start_OutOfFrame",
                                    "De_novo_Start_InFrame", "Missense", "Read-through", "Indel")
      genesWithProteinChangeList <- list()

      #Fill those lists.
      #_______________________________________
      for(i in 1:length(unique)){
        print(paste(i, "of", length(unique)), sep='')
        indices<-unique[i] == sampIds

        #Now for each sample, pull out a list of genes with somatic mutations
        variantType <- data[indices,"Variant_Classification"]
        theGenes <- data[indices,"Hugo_Symbol"]
        names(variantType) <- theGenes
        mutsListAll[[i]] <- variantType #A list of all the mutation occurring in each sample.
        genesWithProteinChangeList[[i]] <- unique(names(variantType[variantType %in% proteinChangingMutations]))
      }

      allMutatedGenes <- unique(names(unlist(mutsListAll))) #All the unique genes mutated across all patients.
      mutationTypes <- table(unlist(mutsListAll)) #Frequency of mutation type across all patients.
      mutsListAll_unlist <- unlist(mutsListAll, recursive=F) #Mutation type and gene association.
      genesWithProteinChangeList_unlist <- unlist(genesWithProteinChangeList, recursive=F) #Genes with protein changing mutations across all patients.

      #From mutsListAll we can then create a matrix indicating if the gene has a coding mutation.
      #_______________________________________
      mutMat <- numeric((length(unique)*length(allMutatedGenes)))
      dim(mutMat) <- c(length(allMutatedGenes), length(unique))
      rownames(mutMat) <- allMutatedGenes
      colnames(mutMat) <- unique

      #Now populate this matrix with the relevant information about what kind of mutation each gene has in each sample.
      #_______________________________________
      for(i in 1:length(unique)){
        print(paste(i, "of", length(unique)), sep='')
        mutMat[genesWithProteinChangeList[[i]], i] <- rep(1, length(genesWithProteinChangeList[[i]]))
      }

      #Identify samples that occur in drug and mutation data.
      #_______________________________________
      drug_samps<-colnames(drug_prediction)
      inPredAndMutData <- drug_samps[drug_samps %in% unique] #Samples for which we have both predicted drug response and mutation calls

      #Run the associations between all genes and drugs, for drugs with at least 50 mutations.
      #_______________________________________
      drug_prediction_filt_ord <- drug_prediction[, inPredAndMutData] #The preds for the 01A samples we have both prediction and mutation data for.
      mutMat_nodups_ordFilt <- mutMat[, inPredAndMutData]
      commonMuts <- apply(mutMat_nodups_ordFilt, 1, sum)
      if (length(which(commonMuts >= n)) == 0){
        stop((paste("\nERROR: Mutations were not identified in at least", n, "genes. Recommend decreasing the n parameter.", sep=" ")))
      }
      commonlyMutated <- mutMat_nodups_ordFilt[which(commonMuts >= n), ]

      #If there are gene entries with an unknown HUGO ID, remove it.
      #_______________________________________
      if("Unknown" %in% rownames(commonlyMutated)){
        indices<-'Unknown' %in% rownames(commonlyMutated)
        commonlyMutated<-commonlyMutated[-indices,]
      }

      #Get p values and beta values.
      #_______________________________________
      pValList <- list()
      betaValList <- list()
      for(i in 1:nrow(drug_prediction_filt_ord)){ #For each drug...
        pValList[[i]] <- numeric()
        betaValList[[i]] <- numeric()
        for(j in 1:nrow(commonlyMutated))
        {
          thecoefs <- coef(summary(lm(drug_prediction_filt_ord[i,]~commonlyMutated[j,])))
          pValList[[i]][[j]] <- thecoefs[2,4]
          betaValList[[i]][[j]] <- thecoefs[2,1]
        }
      }

      #Get the adjusted p-value for each gene-drug combination, pull out the significant associations
      #and create a supplementary table that lists these for "predictable" drugs?.
      #_______________________________________
      sigPs <- list()
      pAdjListCantype <- list()
      for(i in 1:length(pValList))
      {
        names(pValList[[i]]) <- rownames(commonlyMutated)
        names(betaValList[[i]]) <- rownames(commonlyMutated)
        padj <- p.adjust(pValList[[i]], method="BH")
        sigPs[[i]] <- padj[padj < 0.05]
        pAdjListCantype[[i]] <- padj
      }
      names(sigPs) <- rownames(drug_prediction_filt_ord)
      names(pValList) <- rownames(drug_prediction_filt_ord)
      names(betaValList) <- rownames(drug_prediction_filt_ord)
      names(pAdjListCantype) <- rownames(drug_prediction_filt_ord)

      output<-sort(unlist(pValList))[1:30]

      write.table(output, file='./MutationTestOutput.txt')

      #Print the top associations
      #return(output)
    }
  }
}
