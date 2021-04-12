# CRUK_CI_Summer_School_2021_ScRnaSeq
Introduction to single-cell RNA-seq analysis at the CRUK Summer School 2021

**Private repo** to update materials used last year and bundled it as the bookdown initiated last year but not finished.

Creating a bookdown can be done in two ways:

* merging all chapter first then rendering: the 'merge-knit' way
* rendering each chapter separately then merging: the 'knit-merge' way

I prefer the 'knit-merge' way as each chapter is rendered separately, each in a different session. I also like to skip unnecessary re-running during development and so use caching. That however sometimes creates issues when a cached chunk uses an object updated elsewhere (another chunk) since the last caching. I ended up in such a situation and decided to use the 'merge-knit' method instead, with no caching, which also offers the possibility to render chapters one at a time. Which is fine, until parameters are added to a chapter that are not listed in the index template which is the only one used with this single-session method (as far as I understand. But I have not yet checked that by editing index.Rmd). This is currently the case here as I added the 'setSuf' parameter to the data set integration chapters. The latter were thus rendered as stand-alone documents, ie not included in the bookdown document. That is the temporary situation on Mon 12 April 2021.

The repo so far aims to gather in a single place the various analyses performed to develop the materials. Not all analyses will be included in the final document used by students during the course. Analyses of the whole data set for example take too long. Once all analyses are compiled, the bookdown will serve as a source to derive another version to share with students, where some chapters are clearly labelled 'exercise', e.g. data set integration, and/or where whole-data-set analyses were excluded or available as annexes.  

Two data sets:

* 'CaronBourque2020': pediatric leukemia, with four sample types, including:
  * pediatric Bone Marrow Mononuclear Cells (PBMMCs)
  * three tumour types: ETV6-RUNX1, HHD, PRE-T and  
* 'HCA': adult BMMCs obtained from the Human Cell Atlas (HCA)
  * (here downsampled from 25000 to 5000 cells per sample)

List of chapters:

* index: index.Rmd
* sequence quality: seqQual.Rmd
* alignmment and feature counting, with cellranger: cellRanger.Rmd
* pre-processing - all cells: preProcAllCells.Rmd
* pre-processing - downsampled sample set, down to 500 cells each: preProc.Rmd
* normalisation - all cells: normalisation_whole.Rmd
* normalisation - downsampled sample set: normalisation_5hcps.Rmd
* normalisation - GSM3872434 sample only: normalisation_GSM3872434_simple.Rmd
* dimensionality reduction for visualisation: dimRedForViz.Rmd
* identification of confounding factor: confounding_caron.Rmd
* feature selection: featSelec.Rmd
* batch effect: batch_GSM3872442.Rmd
* dimensionality reduction for analysis: dimRedForAna.Rmd
* data sets integration - PBMMC: dataSetIntegration_PBMMC.Rmd
* data sets integration - ETV6-RUNX1: dataSetIntegration_ETV6-RUNX1.Rmd
* data sets integration - HHD: dataSetIntegration_HHD.Rmd
* data sets integration - PRE-T: dataSetIntegration_PRE-T.Rmd
* data sets integration - PBMMC + ETV6-RUNX1: dataSetIntegration_PBMMC_ETV6-RUNX1.Rmd
* clustering - merge data sets: clusteringPostDsi.Rmd
* cluster marker gene identification: clusterMarkerGenes.Rmd
* cell cycle assignment (last year's)
* differential expression between conditions (last year's)
* trajectory analysis (last year's)
* doublet detection (last year's)

TODO HERE:

* add bookdown and link to it
* add code to render bookdown
