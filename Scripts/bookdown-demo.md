--- 
title: "CRUK Bioinformatics Summer School 2020 - single-cell RNA-seq analysis"
subtitle: 'Heterogeneity in childhood acute lymphoblastic leukemia with droplet-based 10X Chromium assay.'
author: "Stephane Ballereau, Zeynep Kalender Atak"
date: "2021-03-04"
site: bookdown::bookdown_site
output:
  bookdown::gitbook:
      toc_depth: 6
      includes:
        in_header: header.html
documentclass: book
bibliography: [book.bib]
biblio-style: apalike
link-citations: yes
github-repo: rstudio/bookdown-demo
description: "A report for all analyses (bookdown::gitbook)."
params:
  outDirBit: "AnaWiSce/Ana1"
---

**WORKING DOCUMENT - IN PROGRESS**

<!--
Useful links: https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/SingleCell/practical/QC_and_normalization.html
-->




```
## Warning: The `path` argument of `write_lines()` is deprecated as of readr 1.4.0.
## Please use the `file` argument instead.
```

# Preamble {#Preamble}



## The study

"Childhood acute lymphoblastic leukemia (cALL) is the most common pediatric cancer. It is characterized by bone marrow lymphoid precursors that acquire genetic alterations, resulting in disrupted maturation and uncontrollable proliferation." [Caron et al. 2020](https://www.nature.com/articles/s41598-020-64929-x). Nowaways, up to 85–90% of patients are cured, but others do not respond to treatment or relapse and die. The aim of the study is to characterise the heterogeneity of gene expression at the cell level, within and between patients.

<!--
See similar findings in CRC
https://www.cell.com/cell-stem-cell/fulltext/S1934-5909(20)30151-X?rss=yes

https://www.nature.com/articles/s41568-020-0276-8
https://www.sciencedirect.com/science/article/pii/S1934590920302034

https://scinapse.io/papers/2955045010
(https://www.nature.com/articles/nm.4505)
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4273416/
https://www.nature.com/articles/s41375-018-0127-8
-->

Four type of samples are considered:

* eight patients:
  * six B-ALL
  * two T-ALL
* three healthy pediatric controls
* eight healthy adult controls, publicly available

As the study aims at identifying cell populations, large numbers of cells were sequenced with the droplet-based 10X Chromium assay.

## The plan

We will follow several steps:

* sequencing quality check
* alignment of reads to the human genome (GRCh38) with 10X software cellranger
* quality control (cell calls, cells and genes filtering) <!-- mention doublet detection -->
* count normalisation
* data set integration
* feature selection
* dimensionality reduction
* clustering
* marker gene identification
* cell type annotation
* cell cycle assignment
* trajectory analysis

<!--
TODO: add analysis flowchart?
-->

## The analyses

This report includes:

* **Sequence Quality** is good (see \@ref(SeqQualTop)) <!-- seqQual.Rmd -->
* **cellranger** output suggests the data sets are good quality (see \@ref(AliFeatCountTop)) <!-- cellRanger.Rmd -->
* **Quality control** - a first glance at the RNA data set comprising all droplets deemed to contain at leat one cell, to check the data quality and biological signal (eg. do we observe cell types we expect?) (see \@ref(PreProcTop)) <!-- preProc.Rmd -->
* an annotation of RNA clusters using known PBMC marker genes (see \@ref(AllDropletsCellTypeAnnotTop)) <!-- .Rmd -->

## Overall summary

<!-- seurat3Adt_glanceAtGex.Rmd -->
* 

<!-- seurat3Adt_10x3kPbmcMarkers.Rmd -->
* 

<!-- seurat3Hto.Rmd -->
* 

## Conclusion



## Abbreviations

* BMMC: Bone Marrow Mononuclear Cell
* HCA: Human Cell Atlas
* PCA: Principal Component Analysis
* UMI: Unique Molecular Identifier





<!--chapter:end:index.Rmd-->

---
title: "CRUK CI Summer School 2020 - introduction to single-cell RNA-seq analysis"
subtitle: 'Sequence quality'

author: "Stephane Ballereau"
output:
  html_notebook:
    code_folding: hide
    toc: yes
    toc_float: yes
    number_sections: true
  html_document:
    df_print: paged
    toc: yes
    number_sections: true
    code_folding: hide
  html_book:
    code_folding: hide
params:
  projDir: "/ssd/personal/baller01/20200511_FernandesM_ME_crukBiSs2020"
  inpDirBit: "AnaWiSce/Ana1"
  outDirBit: "AnaWiSce/Ana1"
---

<!--
TODO:
-->




```r
#projDirOsx <- "/Users/baller01/MyMount/clust1a/20200511_FernandesM_ME_crukBiSs2020"
#projDir <- "/mnt/scratcha/bioinformatics/baller01/20200511_FernandesM_ME_crukBiSs2020"
#projDir <- "/home/ubuntu/Course_Materials/scRNAseq"
projDir <- params$projDir
projDirLink <- "/Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020"
inpDirBit <- params$inpDirBit
outDirBit <- params$outDirBit
```


```r
library(DT)
```

# Sequence Quality {#SeqQualTop}

## Introduction

We will use two sets of Bone Marrow Mononuclear Cells (BMMC):

* 'CaronBourque2020': pediatric samples
* 'Hca': HCA Census of Immune Cells for adult BMMCs

Fastq files were retrieved from publicly available archive (SRA and HCA). 

Sequencing quality was assessed and visualised using fastQC and MultiQC.

Library structure reminder:

<!--
![](/ssd/personal/baller01/20200511_FernandesM_ME_crukBiSs2020/Images/tenxLibStructureV3.png)
-->

* The **sample index** identifies the library, with one I7 index per sample
* The 10X **cell barcode** (or cell index) identifies the droplet in the library
* The **UMI** identifies the transcript molecule within a cell and gene
* The **insert** is the transcript molecule, ie the cDNA sequence

Each sample is described with three sets of fastq files:

* **I1**: sample index
* **R1**: 10x barcode + UMI
* **R2**: insert sequence

The sample index is actually a set of four 8-ntd oligo.
For example SIGAB8 is 'AAAGTGCT-GCTACCTG-TGCTGTAA-CTGCAAGC'.
All four are used and identified by a digit, eg 1-4.
Depending on the processing pipeline, fastq files may be returned for each 8-ntd index, or combined into a single file.

For the Caron data set they are combined in a single file, and files for separate lanes were also combined into a single fastq file.

Each sample is identified by three fastq files, one per read type:

* **sample** _ S0 _ L001 _ **I1** _ 001 _ .fastq.gz: contains sample index
* **sample** _ S0 _ L001 _ **R1** _ 001 _ .fastq.gz: contains 10x barcode + UMI
* **sample** _ S0 _ L001 _ **R2** _ 001 _ .fastq.gz: contains insert sequence

We kept the same names for the fastqc output. With for example sample 'SRR9264343':

* **SRR9264343** _ S0 _ L001 _ **I1** _ 001 _ fastqc.html
* **SRR9264343** _ S0 _ L001 _ **R1** _ 001 _ fastqc.html
* **SRR9264343** _ S0 _ L001 _ **R2** _ 001 _ fastqc.html


```r
fastqcDir <- sprintf("%s/Data/%s/fastqc", projDir, "CaronBourque2020")
fastqcDirLink <- sprintf("%s/Data/%s/fastqc", projDirLink, "CaronBourque2020")
```

## CaronBourque2020 - fastqc


```r
# CaronBourque2020
cb_sampleSheetFn <- file.path(projDir, "Data/CaronBourque2020/SraRunTable.txt")
cb_sampleSheet <- read.table(cb_sampleSheetFn, header=T, sep=",")
#cb_sampleSheet <-  cb_sampleSheet %>% filter(!Run == "SRR9264351")
cb_sampleSheet
```

```
##           Run Assay.Type AvgSpotLen       Bases  BioProject    BioSample
## 1  SRR9264343    RNA-Seq        132 27850288884 PRJNA548203 SAMN12011162
## 2  SRR9264344    RNA-Seq        132 43613421192 PRJNA548203 SAMN12011172
## 3  SRR9264345    RNA-Seq        132 43838527392 PRJNA548203 SAMN12011171
## 4  SRR9264346    RNA-Seq        132 39752529300 PRJNA548203 SAMN12011170
## 5  SRR9264347    RNA-Seq        132 41035092252 PRJNA548203 SAMN12011169
## 6  SRR9264348    RNA-Seq        132 42840756288 PRJNA548203 SAMN12011168
## 7  SRR9264349    RNA-Seq        132 42953865372 PRJNA548203 SAMN12011167
## 8  SRR9264350    RNA-Seq        132 42822420960 PRJNA548203 SAMN12011166
## 9  SRR9264351    RNA-Seq        132 28322630028 PRJNA548203 SAMN12011165
## 10 SRR9264352    RNA-Seq        132 36199482528 PRJNA548203 SAMN12011165
## 11 SRR9264353    RNA-Seq        132 41446760124 PRJNA548203 SAMN12011164
## 12 SRR9264354    RNA-Seq        132 42802129128 PRJNA548203 SAMN12011163
##          Bytes
## 1  18644549905
## 2  27638885644
## 3  28054431102
## 4  25564104997
## 5  24777477094
## 6  27432674292
## 7  27523442193
## 8  27282064655
## 9  19040444664
## 10 22143300246
## 11 26850120365
## 12 27774281557
##                                                            Cell_type
## 1     Pre-B t(12;21) [ETV6-RUNX1] acute lymphoblastic leukemia cells
## 2     Pre-B t(12;21) [ETV6-RUNX1] acute lymphoblastic leukemia cells
## 3     Pre-B t(12;21) [ETV6-RUNX1] acute lymphoblastic leukemia cells
## 4     Pre-B t(12;21) [ETV6-RUNX1] acute lymphoblastic leukemia cells
## 5  Pre-B High hyper diploid [HHD] acute lymphoblastic leukemia cells
## 6  Pre-B High hyper diploid [HHD] acute lymphoblastic leukemia cells
## 7                           Pre-T acute lymphoblastic leukemia cells
## 8                           Pre-T acute lymphoblastic leukemia cells
## 9                    Healthy pediatric bone marrow mononuclear cells
## 10                   Healthy pediatric bone marrow mononuclear cells
## 11                   Healthy pediatric bone marrow mononuclear cells
## 12                   Healthy pediatric bone marrow mononuclear cells
##    Center.Name Consent DATASTORE.filetype DATASTORE.provider
## 1          GEO  public          fastq,sra         gs,ncbi,s3
## 2          GEO  public          fastq,sra         gs,ncbi,s3
## 3          GEO  public          fastq,sra         gs,ncbi,s3
## 4          GEO  public          fastq,sra         gs,ncbi,s3
## 5          GEO  public          fastq,sra         gs,ncbi,s3
## 6          GEO  public          fastq,sra         gs,ncbi,s3
## 7          GEO  public          fastq,sra         gs,ncbi,s3
## 8          GEO  public          fastq,sra         gs,ncbi,s3
## 9          GEO  public          fastq,sra         gs,ncbi,s3
## 10         GEO  public          fastq,sra         gs,ncbi,s3
## 11         GEO  public          fastq,sra         gs,ncbi,s3
## 12         GEO  public          fastq,sra         gs,ncbi,s3
##                  DATASTORE.region                          disease_state
## 1  gs.US,ncbi.public,s3.us-east-1 Childhood acute lymphoblastic leukemia
## 2  gs.US,ncbi.public,s3.us-east-1 Childhood acute lymphoblastic leukemia
## 3  gs.US,ncbi.public,s3.us-east-1 Childhood acute lymphoblastic leukemia
## 4  gs.US,ncbi.public,s3.us-east-1 Childhood acute lymphoblastic leukemia
## 5  gs.US,ncbi.public,s3.us-east-1 Childhood acute lymphoblastic leukemia
## 6  gs.US,ncbi.public,s3.us-east-1 Childhood acute lymphoblastic leukemia
## 7  gs.US,ncbi.public,s3.us-east-1 Childhood acute lymphoblastic leukemia
## 8  gs.US,ncbi.public,s3.us-east-1 Childhood acute lymphoblastic leukemia
## 9  gs.US,ncbi.public,s3.us-east-1              Healthy pediatric control
## 10 gs.US,ncbi.public,s3.us-east-1              Healthy pediatric control
## 11 gs.US,ncbi.public,s3.us-east-1              Healthy pediatric control
## 12 gs.US,ncbi.public,s3.us-east-1              Healthy pediatric control
##    Experiment GEO_Accession..exp.          Instrument LibraryLayout
## 1  SRX6034681          GSM3872434 Illumina HiSeq 4000        PAIRED
## 2  SRX6034682          GSM3872435 Illumina HiSeq 4000        PAIRED
## 3  SRX6034683          GSM3872436 Illumina HiSeq 4000        PAIRED
## 4  SRX6034684          GSM3872437 Illumina HiSeq 4000        PAIRED
## 5  SRX6034685          GSM3872438 Illumina HiSeq 4000        PAIRED
## 6  SRX6034686          GSM3872439 Illumina HiSeq 4000        PAIRED
## 7  SRX6034687          GSM3872440 Illumina HiSeq 4000        PAIRED
## 8  SRX6034688          GSM3872441 Illumina HiSeq 4000        PAIRED
## 9  SRX6034689          GSM3872442 Illumina HiSeq 4000        PAIRED
## 10 SRX6034689          GSM3872442 Illumina HiSeq 4000        PAIRED
## 11 SRX6034690          GSM3872443 Illumina HiSeq 4000        PAIRED
## 12 SRX6034691          GSM3872444 Illumina HiSeq 4000        PAIRED
##    LibrarySelection  LibrarySource     Organism Platform          ReleaseDate
## 1              cDNA TRANSCRIPTOMIC Homo sapiens ILLUMINA 2020-02-14T00:00:00Z
## 2              cDNA TRANSCRIPTOMIC Homo sapiens ILLUMINA 2020-02-14T00:00:00Z
## 3              cDNA TRANSCRIPTOMIC Homo sapiens ILLUMINA 2020-02-14T00:00:00Z
## 4              cDNA TRANSCRIPTOMIC Homo sapiens ILLUMINA 2020-02-14T00:00:00Z
## 5              cDNA TRANSCRIPTOMIC Homo sapiens ILLUMINA 2020-02-14T00:00:00Z
## 6              cDNA TRANSCRIPTOMIC Homo sapiens ILLUMINA 2020-02-14T00:00:00Z
## 7              cDNA TRANSCRIPTOMIC Homo sapiens ILLUMINA 2020-02-14T00:00:00Z
## 8              cDNA TRANSCRIPTOMIC Homo sapiens ILLUMINA 2020-02-14T00:00:00Z
## 9              cDNA TRANSCRIPTOMIC Homo sapiens ILLUMINA 2020-02-14T00:00:00Z
## 10             cDNA TRANSCRIPTOMIC Homo sapiens ILLUMINA 2020-02-14T00:00:00Z
## 11             cDNA TRANSCRIPTOMIC Homo sapiens ILLUMINA 2020-02-14T00:00:00Z
## 12             cDNA TRANSCRIPTOMIC Homo sapiens ILLUMINA 2020-02-14T00:00:00Z
##    Sample.Name source_name SRA.Study
## 1   GSM3872434  ETV6-RUNX1 SRP201012
## 2   GSM3872435  ETV6-RUNX1 SRP201012
## 3   GSM3872436  ETV6-RUNX1 SRP201012
## 4   GSM3872437  ETV6-RUNX1 SRP201012
## 5   GSM3872438         HHD SRP201012
## 6   GSM3872439         HHD SRP201012
## 7   GSM3872440       PRE-T SRP201012
## 8   GSM3872441       PRE-T SRP201012
## 9   GSM3872442       PBMMC SRP201012
## 10  GSM3872442       PBMMC SRP201012
## 11  GSM3872443       PBMMC SRP201012
## 12  GSM3872444       PBMMC SRP201012
```


```r
htmlVec <- list.files(fastqcDir)
htmlVec <- grep("\\.html$", htmlVec, value=TRUE)
```


```r
filesDf <- data.frame(
		      "I1" = sprintf("%s_S0_L001_%s_001_fastqc.html", cb_sampleSheet$Run, "I1"),
		      "R1" = sprintf("%s_S0_L001_%s_001_fastqc.html", cb_sampleSheet$Run, "R1"),
		      "R2" = sprintf("%s_S0_L001_%s_001_fastqc.html", cb_sampleSheet$Run, "R2")
)
rownames(filesDf) <- cb_sampleSheet$Run
```


```r
for (runx in cb_sampleSheet$Run)
{
	cat("Run ", runx, ":\n\n")
	for(i in c("I1", "R1", "R2"))
	{
		#filepath <- file.path(fastqcDir, filesDf[runx,i])
		filepath <- file.path(fastqcDirLink, filesDf[runx,i])
		cat(i, ": [", filesDf[runx,i], "](",filepath,")\n\n")
	}
}
```

Run  SRR9264343 :

I1 : [ SRR9264343_S0_L001_I1_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264343_S0_L001_I1_001_fastqc.html )

R1 : [ SRR9264343_S0_L001_R1_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264343_S0_L001_R1_001_fastqc.html )

R2 : [ SRR9264343_S0_L001_R2_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264343_S0_L001_R2_001_fastqc.html )

Run  SRR9264344 :

I1 : [ SRR9264344_S0_L001_I1_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264344_S0_L001_I1_001_fastqc.html )

R1 : [ SRR9264344_S0_L001_R1_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264344_S0_L001_R1_001_fastqc.html )

R2 : [ SRR9264344_S0_L001_R2_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264344_S0_L001_R2_001_fastqc.html )

Run  SRR9264345 :

I1 : [ SRR9264345_S0_L001_I1_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264345_S0_L001_I1_001_fastqc.html )

R1 : [ SRR9264345_S0_L001_R1_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264345_S0_L001_R1_001_fastqc.html )

R2 : [ SRR9264345_S0_L001_R2_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264345_S0_L001_R2_001_fastqc.html )

Run  SRR9264346 :

I1 : [ SRR9264346_S0_L001_I1_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264346_S0_L001_I1_001_fastqc.html )

R1 : [ SRR9264346_S0_L001_R1_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264346_S0_L001_R1_001_fastqc.html )

R2 : [ SRR9264346_S0_L001_R2_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264346_S0_L001_R2_001_fastqc.html )

Run  SRR9264347 :

I1 : [ SRR9264347_S0_L001_I1_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264347_S0_L001_I1_001_fastqc.html )

R1 : [ SRR9264347_S0_L001_R1_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264347_S0_L001_R1_001_fastqc.html )

R2 : [ SRR9264347_S0_L001_R2_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264347_S0_L001_R2_001_fastqc.html )

Run  SRR9264348 :

I1 : [ SRR9264348_S0_L001_I1_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264348_S0_L001_I1_001_fastqc.html )

R1 : [ SRR9264348_S0_L001_R1_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264348_S0_L001_R1_001_fastqc.html )

R2 : [ SRR9264348_S0_L001_R2_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264348_S0_L001_R2_001_fastqc.html )

Run  SRR9264349 :

I1 : [ SRR9264349_S0_L001_I1_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264349_S0_L001_I1_001_fastqc.html )

R1 : [ SRR9264349_S0_L001_R1_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264349_S0_L001_R1_001_fastqc.html )

R2 : [ SRR9264349_S0_L001_R2_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264349_S0_L001_R2_001_fastqc.html )

Run  SRR9264350 :

I1 : [ SRR9264350_S0_L001_I1_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264350_S0_L001_I1_001_fastqc.html )

R1 : [ SRR9264350_S0_L001_R1_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264350_S0_L001_R1_001_fastqc.html )

R2 : [ SRR9264350_S0_L001_R2_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264350_S0_L001_R2_001_fastqc.html )

Run  SRR9264351 :

I1 : [ SRR9264351_S0_L001_I1_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264351_S0_L001_I1_001_fastqc.html )

R1 : [ SRR9264351_S0_L001_R1_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264351_S0_L001_R1_001_fastqc.html )

R2 : [ SRR9264351_S0_L001_R2_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264351_S0_L001_R2_001_fastqc.html )

Run  SRR9264352 :

I1 : [ SRR9264352_S0_L001_I1_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264352_S0_L001_I1_001_fastqc.html )

R1 : [ SRR9264352_S0_L001_R1_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264352_S0_L001_R1_001_fastqc.html )

R2 : [ SRR9264352_S0_L001_R2_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264352_S0_L001_R2_001_fastqc.html )

Run  SRR9264353 :

I1 : [ SRR9264353_S0_L001_I1_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264353_S0_L001_I1_001_fastqc.html )

R1 : [ SRR9264353_S0_L001_R1_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264353_S0_L001_R1_001_fastqc.html )

R2 : [ SRR9264353_S0_L001_R2_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264353_S0_L001_R2_001_fastqc.html )

Run  SRR9264354 :

I1 : [ SRR9264354_S0_L001_I1_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264354_S0_L001_I1_001_fastqc.html )

R1 : [ SRR9264354_S0_L001_R1_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264354_S0_L001_R1_001_fastqc.html )

R2 : [ SRR9264354_S0_L001_R2_001_fastqc.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc/SRR9264354_S0_L001_R2_001_fastqc.html )

## CaronBourque2020 - MultiQC

### sample index: I1


```r
htmlVec <- list.files(paste0(fastqcDir, "/Multiqc/I1"))
htmlVec <- grep("\\.html$", htmlVec, value=TRUE)
for(i in htmlVec){
	filename <- file.path(fastqcDirLink, "/Multiqc/I1", i)
	cat("[", i, "](",filename,")\n\n")
}
```

[ multiqc_report.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc//Multiqc/I1/multiqc_report.html )

### cell barcode + UMI: R1


```r
htmlVec <- list.files(paste0(fastqcDir, "/Multiqc/R1"))
htmlVec <- grep("\\.html$", htmlVec, value=TRUE)
for(i in htmlVec){
	filename <- file.path(fastqcDirLink, "/Multiqc/R1", i)
	cat("[", i, "](",filename,")\n\n")
}
```

[ multiqc_report.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc//Multiqc/R1/multiqc_report.html )

### insert: R2


```r
htmlVec <- list.files(paste0(fastqcDir, "/Multiqc/R2"))
htmlVec <- grep("\\.html$", htmlVec, value=TRUE)
for(i in htmlVec){
	filename <- file.path(fastqcDirLink, "/Multiqc/R2", i)
	cat("[", i, "](",filename,")\n\n")
}
```

[ multiqc_report.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/CaronBourque2020/fastqc//Multiqc/R2/multiqc_report.html )

## HCA adult BMMC - fastqc

For the HCA adult BMMC fastq files were provided for each 8-ntd sample index and lane. We ran fastqc on each separately. We are therefore not listing links to the fastqc reports but only to the MultiQC reports.


```r
fastqcDir <- sprintf("%s/Data/%s/fastqc", projDir, "Hca")
fastqcDirLink <- sprintf("%s/Data/%s/fastqc", projDirLink, "Hca")

# HCA
hca_sampleSheetFn <- file.path(projDir, "Data/Hca/accList_Hca.txt")

hca_sampleSheet <- read.table(hca_sampleSheetFn, header=F, sep=",")
colnames(hca_sampleSheet) <- "Run"
hca_sampleSheet
```

```
##         Run
## 1 MantonBM1
## 2 MantonBM2
## 3 MantonBM3
## 4 MantonBM4
## 5 MantonBM5
## 6 MantonBM6
## 7 MantonBM7
## 8 MantonBM8
```


```r
htmlVec <- list.files(fastqcDir)
htmlVec <- grep("\\.html$", htmlVec, value=TRUE)
```


378 fastqc reports were compiled in the multiQC reports below.

##  HCA adult BMMC - MultiQC

### sample index: I1


```r
htmlVec <- list.files(paste0(fastqcDir, "/Multiqc/I1"))
htmlVec <- grep("\\.html$", htmlVec, value=TRUE)
for(i in htmlVec){
	filename <- file.path(fastqcDirLink, "/Multiqc/I1", i)
	cat("[", i, "](",filename,")\n\n")
}
```

[ multiqc_report.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/Hca/fastqc//Multiqc/I1/multiqc_report.html )

### cell barcode + UMI: R1


```r
htmlVec <- list.files(paste0(fastqcDir, "/Multiqc/R1"))
htmlVec <- grep("\\.html$", htmlVec, value=TRUE)
for(i in htmlVec){
	filename <- file.path(fastqcDirLink, "/Multiqc/R1", i)
	cat("[", i, "](",filename,")\n\n")
}
```

[ multiqc_report.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/Hca/fastqc//Multiqc/R1/multiqc_report.html )

### insert: R2


```r
htmlVec <- list.files(paste0(fastqcDir, "/Multiqc/R2"))
htmlVec <- grep("\\.html$", htmlVec, value=TRUE)
for(i in htmlVec){
	filename <- file.path(fastqcDirLink, "/Multiqc/R2", i)
	cat("[", i, "](",filename,")\n\n")
}
```

[ multiqc_report.html ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Data/Hca/fastqc//Multiqc/R2/multiqc_report.html )



<!--chapter:end:seqQual.Rmd-->

---
title: "CRUK CI Summer School 2020 - introduction to single-cell RNA-seq analysis"
subtitle: 'cellranger'

author: "Stephane Ballereau"
output:
  html_notebook:
    code_folding: hide
    toc: yes
    toc_float: yes
    number_sections: true
  html_document:
    df_print: paged
    toc: yes
    number_sections: true
    code_folding: hide
  html_book:
    code_folding: hide
params:
  projDir: "/ssd/personal/baller01/20200511_FernandesM_ME_crukBiSs2020"
  inpDirBit: "AnaWiSce/Ana1"
  outDirBit: "AnaWiSce/Ana1"
---

**WORKING DOCUMENT - IN PROGRESS**

<!--
TODO:
-->





```r
library(dplyr)
```

# Alignment and feature counting {#AliFeatCountTop}

## Introduction

We will use two sets of Bone Marrow Mononuclear Cells (BMMC):

* 'CaronBourque2020': pediatric samples
* 'Hca': HCA Census of Immune Cells for adult BMMCs

Fastq files were retrieved from publicly available archive (SRA and HCA). 

Sequencing quality was assessed and visualised using fastQC and MultiQC.

Reads were aligned against GRCh38 and features counted using cellranger (v3.1.0).


```r
#wrkDir <- "/mnt/scratchb/bioinformatics/baller01/20200511_FernandesM_ME_crukBiSs2020/CaronBourque2020/grch38300"
#setwd(wrkDir)
projDir <- params$projDir
#projDirLink <- "/Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020"
projDirLink <- gsub("/ssd/personal/baller01", "/Users/baller01/MyMount/svr008ssd", projDir)
inpDirBit <- params$inpDirBit # "AnaWiSeurat/Attempt1"
outDirBit <- params$outDirBit # "AnaWiSeurat/Attempt1"
plotDir <- "QcPlots"

# eg # cellrangerDirLink <- sprintf("%s/%s/grch38300", projDirLink, "CaronBourque2020")
```

## 10X cellranger pipeline in brief

Each sample was analysed separately with [cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger).
This pipeline "is a set of analysis pipelines that process Chromium single-cell RNA-seq output to align reads, generate feature-barcode matrices and perform clustering and gene expression analysis."

**TODO** Add code to call cellranger

## sample sheet


```r
# CaronBourque2020
cb_sampleSheetFn <- file.path(projDir, "Data/CaronBourque2020/SraRunTable.txt")
# Human Cell Atlas
hca_sampleSheetFn <- file.path(projDir, "Data/Hca/accList_Hca.txt")

# read sample sheet in:
splShtColToKeep <- c("Run", "Sample.Name", "source_name")

cb_sampleSheet <- read.table(cb_sampleSheetFn, header=T, sep=",")
hca_sampleSheet <- read.table(hca_sampleSheetFn, header=F, sep=",")
colnames(hca_sampleSheet) <- "Sample.Name"
hca_sampleSheet$Run <- hca_sampleSheet$Sample.Name
hca_sampleSheet$source_name <- "ABMMC" # adult BMMC

sampleSheetCat <- rbind(cb_sampleSheet[,splShtColToKeep], hca_sampleSheet[,splShtColToKeep])
```


```r
sampleSheetCat %>%
	#DT::datatable(options = list(dom='t'))
	DT::datatable(options = list(pageLength = 10))
```

```{=html}
<div id="htmlwidget-09f9d204e28ebef00012" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-09f9d204e28ebef00012">{"x":{"filter":"none","data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"],["SRR9264343","SRR9264344","SRR9264345","SRR9264346","SRR9264347","SRR9264348","SRR9264349","SRR9264350","SRR9264351","SRR9264352","SRR9264353","SRR9264354","MantonBM1","MantonBM2","MantonBM3","MantonBM4","MantonBM5","MantonBM6","MantonBM7","MantonBM8"],["GSM3872434","GSM3872435","GSM3872436","GSM3872437","GSM3872438","GSM3872439","GSM3872440","GSM3872441","GSM3872442","GSM3872442","GSM3872443","GSM3872444","MantonBM1","MantonBM2","MantonBM3","MantonBM4","MantonBM5","MantonBM6","MantonBM7","MantonBM8"],["ETV6-RUNX1","ETV6-RUNX1","ETV6-RUNX1","ETV6-RUNX1","HHD","HHD","PRE-T","PRE-T","PBMMC","PBMMC","PBMMC","PBMMC","ABMMC","ABMMC","ABMMC","ABMMC","ABMMC","ABMMC","ABMMC","ABMMC"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Run<\/th>\n      <th>Sample.Name<\/th>\n      <th>source_name<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"pageLength":10,"order":[],"autoWidth":false,"orderClasses":false,"columnDefs":[{"orderable":false,"targets":0}]}},"evals":[],"jsHooks":[]}</script>
```

## 10X cellranger reports for CaronBourque2020


```r
#cellrangerDir <- sprintf("%s/%s/grch38300", projDir, "CaronBourque2020")
#projDirOsx <- "/Users/baller01/MyMount/clust1b/20200511_FernandesM_ME_crukBiSs2020"

# make dir name for each sample of interest
# with 'Run' column

sampleSheet <- sampleSheetCat %>%
	filter(! source_name == "ABMMC")

cellrangerDirLink <- sprintf("%s/%s/grch38300", projDirLink, "CaronBourque2020")
htmlVec <- sprintf("%s/%s/%s/outs/web_summary.html", cellrangerDirLink, sampleSheet$Run, sampleSheet$Run)
names(htmlVec) <- sampleSheet$Run

for(i in 1:length(htmlVec)){
	cat("[", names(htmlVec)[i], "](", htmlVec[i],")\n\n")
}
```

[ SRR9264343 ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/CaronBourque2020/grch38300/SRR9264343/SRR9264343/outs/web_summary.html )

[ SRR9264344 ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/CaronBourque2020/grch38300/SRR9264344/SRR9264344/outs/web_summary.html )

[ SRR9264345 ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/CaronBourque2020/grch38300/SRR9264345/SRR9264345/outs/web_summary.html )

[ SRR9264346 ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/CaronBourque2020/grch38300/SRR9264346/SRR9264346/outs/web_summary.html )

[ SRR9264347 ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/CaronBourque2020/grch38300/SRR9264347/SRR9264347/outs/web_summary.html )

[ SRR9264348 ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/CaronBourque2020/grch38300/SRR9264348/SRR9264348/outs/web_summary.html )

[ SRR9264349 ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/CaronBourque2020/grch38300/SRR9264349/SRR9264349/outs/web_summary.html )

[ SRR9264350 ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/CaronBourque2020/grch38300/SRR9264350/SRR9264350/outs/web_summary.html )

[ SRR9264351 ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/CaronBourque2020/grch38300/SRR9264351/SRR9264351/outs/web_summary.html )

[ SRR9264352 ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/CaronBourque2020/grch38300/SRR9264352/SRR9264352/outs/web_summary.html )

[ SRR9264353 ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/CaronBourque2020/grch38300/SRR9264353/SRR9264353/outs/web_summary.html )

[ SRR9264354 ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/CaronBourque2020/grch38300/SRR9264354/SRR9264354/outs/web_summary.html )

```r
# TODO: add links to sample sheet and show with DT::datatable
```

## 10X cellranger reports for HCA's adult BMMCs


```r
#cellrangerDir <- sprintf("%s/%s/grch38300", projDir, "CaronBourque2020")
#projDirOsx <- "/Users/baller01/MyMount/clust1b/20200511_FernandesM_ME_crukBiSs2020"

# make dir name for each sample of interest
# with 'Run' column

sampleSheet <- sampleSheetCat %>%
	filter(source_name == "ABMMC")

cellrangerDirLink <- sprintf("%s/%s/grch38300", projDirLink, "Hca")
htmlVec <- sprintf("%s/%s/%s/outs/web_summary.html", cellrangerDirLink, sampleSheet$Run, sampleSheet$Run)
names(htmlVec) <- sampleSheet$Run

for(i in 1:length(htmlVec)){
	cat("[", names(htmlVec)[i], "](", htmlVec[i],")\n\n")
}
```

[ MantonBM1 ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Hca/grch38300/MantonBM1/MantonBM1/outs/web_summary.html )

[ MantonBM2 ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Hca/grch38300/MantonBM2/MantonBM2/outs/web_summary.html )

[ MantonBM3 ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Hca/grch38300/MantonBM3/MantonBM3/outs/web_summary.html )

[ MantonBM4 ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Hca/grch38300/MantonBM4/MantonBM4/outs/web_summary.html )

[ MantonBM5 ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Hca/grch38300/MantonBM5/MantonBM5/outs/web_summary.html )

[ MantonBM6 ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Hca/grch38300/MantonBM6/MantonBM6/outs/web_summary.html )

[ MantonBM7 ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Hca/grch38300/MantonBM7/MantonBM7/outs/web_summary.html )

[ MantonBM8 ]( /Users/baller01/MyMount/svr008ssd/20200511_FernandesM_ME_crukBiSs2020/Hca/grch38300/MantonBM8/MantonBM8/outs/web_summary.html )

<!--chapter:end:cellRanger.Rmd-->

