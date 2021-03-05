#'This function maps cnv data to genes.
#'The output of this function is a .RData file called map.RData; this file contains theCnvQuantVecList_mat (rows are genes, and columns are samples) and tumorSamps (indicates which samples are primary tumor samples, 01A).
#'@param theRootDir The directory/location in which you downloaded the data to.The output of this function will be saved to this directory.
#'@param Cnvs The output obtained from the download() function (the cnv data); a table with colnames() Sample (named using the TCGA patient barcode), Chromosome, Start, End, Num_Probes, and Segment_Mean.
#'@keywords Map CNV data to genes
#'@return A .RData file called, map.RData, which stores two objects: theCnvQuantVecList_mat (rows are genes, columns are samples), tumorSamps (indicates which samples are primary tumor/01A).
#'@import org.Hs.eg.db
#'@import GenomicFeatures
#'@export
map_cnv<-function(theRootDir, Cnvs){

  #Check colnames() of the cnv data.
  #This is important because depending on what method you use to obtain the cnv data, colnames() might slightly differ (aka Segment.Mean and not Segment_Mean).
  if (('Segment_Mean' %in% colnames(Cnvs)) == FALSE)
    stop("\nERROR: Check colnames() of cnv data. colnames() must include Sample, Chromosome, Start, End, and Segment_Mean")

  #Load the gene ranges for HG19 using.
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  geneRanges <- genes(txdb) #seqnames (chr), ranges, gene id.
  e2s = toTable(org.Hs.egSYMBOL) #gene id, gene symbol.
  syms <- e2s[, "symbol"] #Gene symbols for HG19.
  names(syms) <- e2s[, "gene_id"] #Gene symbols with gene id.
  theGeneSymsOrd <- syms[as.character(geneRanges$gene_id)]

  #path<-paste(theRootDir, '/mut.txt', sep="")
  #Cnvs<-read.table(path, header=TRUE, row.names=1) #The cnv data downloaded with columns sample, chromosome, start, end, num probes, segment mean.
  CnvsList <- split(Cnvs, Cnvs[, "Sample"]) #Put each patient data into list such that each element in the list is for a patient.

  cnvEdGenesList <- list()
  ampGenesList <- list()
  delGenesList <- list()
  numGenesQuantifid <- numeric()
  theCnvQuantVecList <- list()

  for(i in 1:length(CnvsList)){ #Iterate for each patient in the list CnvsList

    chrs <- paste("chr", CnvsList[[i]]$Chromosome, sep="") #chrs contains all the chromosome numbers for these samples in patient i.
    starts <- CnvsList[[i]]$Start #starts contains all the starts.
    ends <- CnvsList[[i]]$End #ends contains all the ends.
    grCnvs <- GRanges(seqnames=Rle(chrs),ranges=IRanges(starts, ends), segMeans=CnvsList[[i]]$Segment_Mean) #chr start-end

    #Amp or del
    grCnvs_ampDel <- grCnvs[grCnvs$segMeans > 1 | grCnvs$segMeans < -1]
    cnvedGenes <- subsetByOverlaps(geneRanges, grCnvs_ampDel, type="within")
    cnvEdGenesList[[i]] <- cnvedGenes$gene_sym

    #Amps
    grCnvs_amp <- grCnvs[grCnvs$segMeans > 1]
    ampedGenes <- subsetByOverlaps(geneRanges, grCnvs_amp, type="within")
    ampGenesList[[i]] <- ampedGenes$gene_sym

    #Dels
    grCnvs_Del <- grCnvs[grCnvs$segMeans < -1]
    deledGenes <- subsetByOverlaps(geneRanges, grCnvs_Del, type="within")
    delGenesList[[i]] <- deledGenes$gene_sym

    #Continuous gene level
    #Use count overlaps to find genes that unambiguously overlap a single peak. Give it an NA it it doesn't overlap a single peak. Assign it the value of the peak if it unambiguously overlaps a peak. PC.
    numOverlaps <- countOverlaps(geneRanges, grCnvs)
    numGenesQuantifid[i] <- sum(numOverlaps == 1)
    inCnv <- which(numOverlaps == 1) # take only gene unambiguously overlaping a peak, this is usually most genes.

    theCnvQuantVec <- rep(NA, length(geneRanges))
    olaps <- findOverlaps(geneRanges, grCnvs, type="within")

    #theCnvQuantVec[olaps@queryHits] <- grCnvs$segMeans[olaps@subjectHits]

    #No slot of name 'subjectHits' for this object of class 'SortedByQueryHits'
    #This problem is caused by a small change in the GenomicAlignments package.
    theCnvQuantVec[olaps@from] <- grCnvs$segMeans[olaps@to]

    theCnvQuantVecList[[i]] <- theCnvQuantVec
    #names(theCnvQuantVecList[[i]]) <- geneRanges$gene_sym
    names(theCnvQuantVecList[[i]]) <- theGeneSymsOrd
  }

  names(theCnvQuantVecList) <- names(CnvsList)
  theCnvQuantVecList_mat <- do.call(rbind, theCnvQuantVecList)
  siteVec <- sapply(strsplit(names(CnvsList), "-"), function(l)return(l[4]))
  tumorSamps <- which(siteVec == "01A")
  save(theCnvQuantVecList_mat, tumorSamps, file=paste(theRootDir, '/', "map.RData", sep="")) # Save these RData files for use by other scripts.

}
