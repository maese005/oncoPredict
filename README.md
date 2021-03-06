# pRRophetic_plus

An R package for drug response prediction and drug-gene association prediction. The prepared GDSC and CTRP matrices for the calcPhenotype() are here: https://drive.google.com/drive/folders/1qgLcsREeGJg_sCYZpTcso2JkWZ08Uu5r?usp=sharing
 *  For drug response prediction, use **calcPhenotype**. This code is based on the paper: pRRophetic: An R Package for Prediction of Clinical Chemotherapeutic Response from Tumor Gene Expression Levels (plos.org)
 *  For pre-clinical biomarker discovery, use **GLDS**. This code is based on the paper: Cancer biomarker discovery is improved by accounting for variability in general levels of drug sensitivity in pre-clinical models | Genome Biology | Full Text (biomedcentral.com
 * For clinical biomarker discovery, use **IDWAS** or indicate **cc=TRUE** in calcPhenotype(). The IDWAS code is based on the paper: Discovering novel pharmacogenomic biomarkers by imputing drug response in cancer patients from large genomics studies (nih.gov)

[ABC](http://example.com)

## R <h2>
 * This directory contains all the R functions included in this package. 

## vignettes <h2> 
  *  This directory contains vignettes which display detailed examples of the functionalities available in this package.
  *  **IDWAS** This directory contains examples of IDWAS code application for clinical drug-gene association prediction. 
      + **cnv.Rmd** Example as to how to download CNV (copy number variation) data from the GDC database, then apply map_cnv() and test().
      + **mut.Rmd** Example as to how to download stomatic mutation data from the GDC database, then apply test(). 

  * **GLDS** This directory contains examples of GLDS code application for pre-clinical drug-gene association prediction. 
      + **glds_GDSC.Rmd** Example of GLDS application to GDSC data. 

  * **calcPhenotype.Rmd** Example of calcPhenotype() application.

## vignetteData <h2>
  * This directory contains the data referenced in the examples provided in the vignette directory. 

## man <h2>
 * This directory contains .Rd (R documentation) files for each function. These files were automatically generated upon creation of the package. 

## NAMESPACE <h2>
 * This file lists the functions to be imported and exported from this package. 

## DESCRIPTION <h2>
 * This file contains the description documentation and metadata for this package. 
