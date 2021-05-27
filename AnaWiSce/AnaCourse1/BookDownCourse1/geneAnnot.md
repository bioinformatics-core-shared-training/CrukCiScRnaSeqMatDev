---
title: "CRUK CI Summer School 2020 - introduction to single-cell RNA-seq analysis"
subtitle: 'Gene Annotation'
author: "Stephane Ballereau"
output:
  html_document:
    df_print: paged
    toc: yes
    number_sections: true
    code_folding: hide
  html_notebook:
    code_folding: hide
    toc: yes
    toc_float: yes
    number_sections: true
  html_book:
    code_folding: hide
params:
  projDir: "/ssd/personal/baller01/20200511_FernandesM_ME_crukBiSs2020"
  inpDirBit: "AnaWiSce/Ana1"
  outDirBit: "AnaWiSce/Ana1"
  bookType: "mk"
  cacheBool: FALSE
---



# Gene annotation {#GeneAnnotTop}

## Introduction

The 10X cellranger pipeline describes gene with Ensembl IDs and gene symbols.
We will need chromosomal location to identify mitochondrial genes to use in quality control.

The online resource used to annotate genes is not always available (see biomaRt package).
We will therefore annotate genes now in a separate chapter that can be run separately.
The annotation is written to a file that is used in following chapters. 

As cellranger generates count matrices with a fixed gene list,
we will use the outcome for a single sample, eg SRR9264343.


```r
projDir <- params$projDir
wrkDir <- sprintf("%s/Data/CaronBourque2020/grch38300", projDir) 
outDirBit <- params$outDirBit
biomartBool <- TRUE # biomaRt sometimes fails, do it once, write to file and use that copy.
```




```r
dir.create(sprintf("%s/%s/Robjects",
		   projDir, outDirBit),
	   showWarnings = FALSE)
```

## Sample sheet

We will load both the Caron and Hca data sets.


```r
# CaronBourque2020
cb_sampleSheetFn <- file.path(projDir, "Data/CaronBourque2020/SraRunTable.txt")
# Human Cell Atlas
hca_sampleSheetFn <- file.path(projDir, "Data/Hca/accList_Hca.txt")

# read sample sheet in:
splShtColToKeep <- c("Run", "Sample.Name", "source_name")

cb_sampleSheet <- read.table(cb_sampleSheetFn, header=T, sep=",")
hca_sampleSheet <- read.table(hca_sampleSheetFn, header=F, sep=",")
colnames(hca_sampleSheet) <- "Sample.Name"
hca_sampleSheet$Run <- hca_sampleSheet$Sample.Name
hca_sampleSheet$source_name <- "ABMMC" # adult BMMC

sampleSheet <- rbind(cb_sampleSheet[,splShtColToKeep], hca_sampleSheet[,splShtColToKeep])

sampleSheet %>%
	as.data.frame() %>%
	DT::datatable(rownames = FALSE)
```

```{=html}
<div id="htmlwidget-8c838d6122beb8c57b8c" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-8c838d6122beb8c57b8c">{"x":{"filter":"none","data":[["SRR9264343","SRR9264344","SRR9264345","SRR9264346","SRR9264347","SRR9264348","SRR9264349","SRR9264350","SRR9264351","SRR9264352","SRR9264353","SRR9264354","MantonBM1","MantonBM2","MantonBM3","MantonBM4","MantonBM5","MantonBM6","MantonBM7","MantonBM8"],["GSM3872434","GSM3872435","GSM3872436","GSM3872437","GSM3872438","GSM3872439","GSM3872440","GSM3872441","GSM3872442","GSM3872442","GSM3872443","GSM3872444","MantonBM1","MantonBM2","MantonBM3","MantonBM4","MantonBM5","MantonBM6","MantonBM7","MantonBM8"],["ETV6-RUNX1","ETV6-RUNX1","ETV6-RUNX1","ETV6-RUNX1","HHD","HHD","PRE-T","PRE-T","PBMMC","PBMMC","PBMMC","PBMMC","ABMMC","ABMMC","ABMMC","ABMMC","ABMMC","ABMMC","ABMMC","ABMMC"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>Run<\/th>\n      <th>Sample.Name<\/th>\n      <th>source_name<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

## SRR9264343

We will load the data for the first sample in the sample sheet: SRR9264343.


```r
i <- 1
sample.path <- sprintf("%s/%s/%s/outs/raw_feature_bc_matrix/",
                       wrkDir,
                       sampleSheet[i,"Run"],
                       sampleSheet[i,"Run"])
sce.raw <- read10xCounts(sample.path, col.names=TRUE)
sce.raw
```

```
## class: SingleCellExperiment 
## dim: 33538 737280 
## metadata(1): Samples
## assays(1): counts
## rownames(33538): ENSG00000243485 ENSG00000237613 ... ENSG00000277475
##   ENSG00000268674
## rowData names(3): ID Symbol Type
## colnames(737280): AAACCTGAGAAACCAT-1 AAACCTGAGAAACCGC-1 ...
##   TTTGTCATCTTTAGTC-1 TTTGTCATCTTTCCTC-1
## colData names(2): Sample Barcode
## reducedDimNames(0):
## altExpNames(0):
```

<!-- Gene annotation -->
<!-- Should be in separate script to run once only -->
<!-- SRR9264343 -->
<!-- raw or filtered: same genes --> 


```r
# retrieve the feature information
gene.info <- rowData(sce.raw)

# setup the biomaRt connection to Ensembl using the correct species genome (hsapiens_gene_ensembl)
ensembl <- useEnsembl(biomart='ensembl',
                      dataset='hsapiens_gene_ensembl',
                      mirror = "www")

#ensembl = useMart(biomart="ensembl",
#		  dataset='hsapiens_gene_ensembl',
#		  host = "www.ensembl.org") # ensemblRedirect = FALSE

# retrieve the attributes of interest from biomaRt using the Ensembl gene ID as the key
# beware that this will only retrieve information for matching IDs
gene_symbol <- getBM(attributes=c('ensembl_gene_id',
                                  'external_gene_name',
                                  'chromosome_name',
                                  'start_position',
                                  'end_position',
                                  'strand'),
                     filters='ensembl_gene_id',
                     mart=ensembl,
                     values=gene.info[, 1])

# create a new data frame of the feature information
gene.merge <- merge(gene_symbol, gene.info, by.x=c('ensembl_gene_id'), by.y=c('ID'), all.y=TRUE)
rownames(gene.merge) <- gene.merge$ensembl_gene_id

# set the order for the same as the original gene information
gene.merge <- gene.merge[gene.info[, 1], ]

# reset the rowdata on the SCE object to contain all of this information
rowData(sce.raw) <- gene.merge
```

Top 10 rows of the gene annotation table:


```r
rowData(sce.raw)[1:10,] %>%
	as.data.frame() %>%
	DT::datatable(rownames = FALSE)
```

```{=html}
<div id="htmlwidget-408ba392df3613304ca1" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-408ba392df3613304ca1">{"x":{"filter":"none","data":[["ENSG00000243485","ENSG00000237613","ENSG00000186092","ENSG00000238009","ENSG00000239945","ENSG00000239906","ENSG00000241599","ENSG00000236601","ENSG00000284733","ENSG00000235146"],["MIR1302-2HG","FAM138A","OR4F5","","","","","","OR4F29",""],["1","1","1","1","1","1","1","1","1","1"],[29554,34554,65419,89295,89551,139790,160446,358857,450740,587629],[31109,36081,71585,133723,91105,140339,161525,366052,451678,594768],[1,-1,1,-1,-1,-1,1,1,-1,1],["MIR1302-2HG","FAM138A","OR4F5","AL627309.1","AL627309.3","AL627309.2","AL627309.4","AL732372.1","OR4F29","AC114498.1"],["Gene Expression","Gene Expression","Gene Expression","Gene Expression","Gene Expression","Gene Expression","Gene Expression","Gene Expression","Gene Expression","Gene Expression"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>ensembl_gene_id<\/th>\n      <th>external_gene_name<\/th>\n      <th>chromosome_name<\/th>\n      <th>start_position<\/th>\n      <th>end_position<\/th>\n      <th>strand<\/th>\n      <th>Symbol<\/th>\n      <th>Type<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```



```r
# slow
# Write object to file
tmpFn <- sprintf("%s/%s/Robjects/sceRaw_postBiomart.Rds", projDir, outDirBit)
saveRDS(sce.raw, tmpFn)
```


```r
# Read object in:
tmpFn <- sprintf("%s/%s/Robjects/sceRaw_postBiomart.Rds", projDir, outDirBit)
sce.raw <- readRDS(tmpFn)
```

## Session information


```r
sessionInfo()
```

```
## R version 4.0.3 (2020-10-10)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: CentOS Linux 8
## 
## Matrix products: default
## BLAS:   /opt/R/R-4.0.3/lib64/R/lib/libRblas.so
## LAPACK: /opt/R/R-4.0.3/lib64/R/lib/libRlapack.so
## 
## locale:
##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] dplyr_1.0.6                 biomaRt_2.46.3             
##  [3] DropletUtils_1.10.3         SingleCellExperiment_1.12.0
##  [5] SummarizedExperiment_1.20.0 Biobase_2.50.0             
##  [7] GenomicRanges_1.42.0        GenomeInfoDb_1.26.7        
##  [9] IRanges_2.24.1              S4Vectors_0.28.1           
## [11] BiocGenerics_0.36.1         MatrixGenerics_1.2.1       
## [13] matrixStats_0.58.0          knitr_1.33                 
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-7              bit64_4.0.5              
##  [3] progress_1.2.2            httr_1.4.2               
##  [5] tools_4.0.3               bslib_0.2.5              
##  [7] utf8_1.2.1                R6_2.5.0                 
##  [9] DT_0.18                   HDF5Array_1.18.1         
## [11] DBI_1.1.1                 rhdf5filters_1.2.0       
## [13] withr_2.4.2               tidyselect_1.1.1         
## [15] prettyunits_1.1.1         bit_4.0.4                
## [17] curl_4.3.1                compiler_4.0.3           
## [19] xml2_1.3.2                DelayedArray_0.16.3      
## [21] bookdown_0.22             sass_0.4.0               
## [23] askpass_1.1               rappdirs_0.3.3           
## [25] stringr_1.4.0             digest_0.6.27            
## [27] rmarkdown_2.8             R.utils_2.10.1           
## [29] XVector_0.30.0            pkgconfig_2.0.3          
## [31] htmltools_0.5.1.1         sparseMatrixStats_1.2.1  
## [33] dbplyr_2.1.1              fastmap_1.1.0            
## [35] limma_3.46.0              htmlwidgets_1.5.3        
## [37] rlang_0.4.11              RSQLite_2.2.7            
## [39] DelayedMatrixStats_1.12.3 jquerylib_0.1.4          
## [41] generics_0.1.0            jsonlite_1.7.2           
## [43] crosstalk_1.1.1           BiocParallel_1.24.1      
## [45] R.oo_1.24.0               RCurl_1.98-1.3           
## [47] magrittr_2.0.1            GenomeInfoDbData_1.2.4   
## [49] scuttle_1.0.4             Matrix_1.3-3             
## [51] Rcpp_1.0.6                Rhdf5lib_1.12.1          
## [53] fansi_0.4.2               lifecycle_1.0.0          
## [55] R.methodsS3_1.8.1         stringi_1.6.1            
## [57] yaml_2.2.1                edgeR_3.32.1             
## [59] zlibbioc_1.36.0           rhdf5_2.34.0             
## [61] BiocFileCache_1.14.0      grid_4.0.3               
## [63] blob_1.2.1                dqrng_0.3.0              
## [65] crayon_1.4.1              lattice_0.20-44          
## [67] beachmat_2.6.4            hms_1.0.0                
## [69] locfit_1.5-9.4            pillar_1.6.1             
## [71] codetools_0.2-18          XML_3.99-0.6             
## [73] glue_1.4.2                evaluate_0.14            
## [75] vctrs_0.3.8               openssl_1.4.4            
## [77] purrr_0.3.4               assertthat_0.2.1         
## [79] cachem_1.0.5              xfun_0.23                
## [81] tibble_3.1.2              AnnotationDbi_1.52.0     
## [83] memoise_2.0.0             ellipsis_0.3.2
```

