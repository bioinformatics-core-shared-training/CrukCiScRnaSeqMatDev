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
  bookType: "mk"
  cacheBool: FALSE
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
<div id="htmlwidget-4a73f7c4e8e51b674319" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-4a73f7c4e8e51b674319">{"x":{"filter":"none","data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"],["SRR9264343","SRR9264344","SRR9264345","SRR9264346","SRR9264347","SRR9264348","SRR9264349","SRR9264350","SRR9264351","SRR9264352","SRR9264353","SRR9264354","MantonBM1","MantonBM2","MantonBM3","MantonBM4","MantonBM5","MantonBM6","MantonBM7","MantonBM8"],["GSM3872434","GSM3872435","GSM3872436","GSM3872437","GSM3872438","GSM3872439","GSM3872440","GSM3872441","GSM3872442","GSM3872442","GSM3872443","GSM3872444","MantonBM1","MantonBM2","MantonBM3","MantonBM4","MantonBM5","MantonBM6","MantonBM7","MantonBM8"],["ETV6-RUNX1","ETV6-RUNX1","ETV6-RUNX1","ETV6-RUNX1","HHD","HHD","PRE-T","PRE-T","PBMMC","PBMMC","PBMMC","PBMMC","ABMMC","ABMMC","ABMMC","ABMMC","ABMMC","ABMMC","ABMMC","ABMMC"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Run<\/th>\n      <th>Sample.Name<\/th>\n      <th>source_name<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"pageLength":10,"order":[],"autoWidth":false,"orderClasses":false,"columnDefs":[{"orderable":false,"targets":0}]}},"evals":[],"jsHooks":[]}</script>
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
