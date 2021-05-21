--- 
title: "CRUK Bioinformatics Summer School 2020 - single-cell RNA-seq analysis"
subtitle: 'Heterogeneity in childhood acute lymphoblastic leukemia with droplet-based 10X Chromium assay.'
author: "Stephane Ballereau, Zeynep Kalender Atak"
date: "2021-05-21"
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
  projDir: "/ssd/personal/baller01/20200511_FernandesM_ME_crukBiSs2020"
  inpDirBit: "AnaWiSce/Ana1"
  outDirBit: "AnaWiSce/Ana1"
  bookType: "mk"
  cacheBool: FALSE
  splSetToGet: "dummy"
  setName: "dummy"
  setSuf: "dummy"
  dsiSuf: 'dummy'
  dirRel: '../..'
---


<!--
  splSetToGet: "PBMMC,ETV6-RUNX1"
  setName: "caron"
  setSuf: "_5hCps"
  dsiSuf: '_dsi'
-->
  
<!--
Useful links: https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/SingleCell/practical/QC_and_normalization.html
-->




```
## Warning: The `path` argument of `write_lines()` is deprecated as of readr 1.4.0.
## Please use the `file` argument instead.
```


<style>
div.blue {background-color:#e6f0ff; border-radius: 5px; padding: 20px;}
</style>
<div class = "blue">

**WORKING DOCUMENT - IN PROGRESS**

</div>



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
    * four 't(12;21)' or 'ETV6-RUNX1'
    * two 'High hyper diploid' or 'HHD'
  * two T-ALL ('PRE-T')
* three healthy pediatric controls
* eight healthy adult controls, publicly available

As the study aims at identifying cell populations, large numbers of cells were sequenced with the droplet-based 10X Chromium assay.

## The plan

We will follow several steps:

* sequencing quality check (see \@ref(SeqQualTop)) <!-- seqQual.Rmd -->
* alignment of reads to the human genome (GRCh38) with 10X software cellranger (see \@ref(AliFeatCountTop)) <!-- cellRanger.Rmd -->
* quality control (cell calls, cells and genes filtering)  (see \@ref(PreProcAllCellsTop) for the 'all-cells' analysis and \@ref(PreProcTop) for the analysis of the downsampled data set to use in the course) <!-- TODO mention doublet detection -->
* count normalisation (see \@ref(NormalisationAllCellsTop) for the 'all-cells' analysis and \@ref(NormalisationCaron5hcpsTop) for the analysis of the downsampled data set to use in the course)
* data set integration (see \@ref(dsi-allCells-PBMMC-ETV6-RUNX1Top) and \@ref(dsi-5hCellPerSpl-PBMMC-ETV6-RUNX1Top), and  \@ref(dsi-allCells-allSetsTop) and \@ref(dsi-5hCellPerSpl-allSetsTop)  )
* feature selection (see \@ref(featSelecTop))
* dimensionality reduction (see \@ref(dimRedForVizTop) for visualisation and \@ref(dimRedForAnaTop) for analysis)
* clustering (see \@ref(clusteringPBMMC-ETV6-RUNX1-5hCellPerSplTop) or \@ref(clustering-with-PBMMC-ETV6-RUNX1-and--5hCellPerSpl))
* marker gene identification (see \@ref(clusterMarkerGenesTop))
* cell cycle assignment (see \@ref(cellCyclePhasesTop))
* cell type annotation (see \@ref(cellTypeAnnotation))
* trajectory analysis (see \@ref(pseudoTimeTop))

<!--
TODO: add analysis flowchart?
-->

## Abbreviations

* ALL: Acute Lymphoblastic Leukemia
* BMMC: Bone Marrow Mononuclear Cell
* cALL: childhood ALL
* HCA: Human Cell Atlas
* PCA: Principal Component Analysis
* UMI: Unique Molecular Identifier

<!-- SETUP -->






