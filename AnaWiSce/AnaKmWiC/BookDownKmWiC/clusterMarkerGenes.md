---
title: "CRUK CI Summer School 2020 - introduction to single-cell RNA-seq analysis"
subtitle: 'Cluster marker genes'

author: "Stephane Ballereau, Zeynep Kalender Atak, Katarzyna Kania"
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
  dirRel: ".."
  inpDirBit: "AnaWiSce/Ana1"
  outDirBit: "AnaWiSce/Ana1"
  bookType: "mk"
  cacheBool: FALSE  
  setName: "caron"
  splSetToGet: "PBMMC,ETV6-RUNX1"
  setSuf: "_5hCellPerSpl"
  dsiSuf: '_dsi'
---



# Cluster marker genes

<img src="../../Images/Andrews2017_Fig1.png" style="margin:auto; display:block" />


```r
if(interactive()) {
  paramsToUse <- params2
} else {
  paramsToUse <- params
}

projDir <- paramsToUse$projDir
dirRel <- paramsToUse$dirRel
outDirBit <- paramsToUse$outDirBit
cacheBool <- paramsToUse$cacheBool
setName <- paramsToUse$setName
splSetToGet <- paramsToUse$splSetToGet
setSuf <- paramsToUse$setSuf
dsiSuf <- paramsToUse$dsiSuf # 'dsi' for data set integration
if(paramsToUse$bookType == "mk"){
	setName <- "caron"
	splSetToGet <- "PBMMC,ETV6-RUNX1"
	setSuf <- "_5hCellPerSpl"
	dsiSuf <- '_dsi'
}
splSetVec <- unlist(strsplit(splSetToGet, ",")) # params may not be read in if knitting book.
splSetToGet2 <- gsub(",", "_", splSetToGet)
nbPcToComp <- 50
figSize <- 7
```




```r
library(ggplot2)
library(scater)
library(scran)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(Cairo)
library(glue)

fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
```

Source: we will follow the [OSCA chapter on marker detection](https://bioconductor.org/books/release/OSCA/marker-detection.html) (with some of its text copied here with little modification). See also the Hemberg group chapter on [differential analysis section](https://scrnaseq-course.cog.sanger.ac.uk/website/biological-analysis.html#dechapter).

To interpret our clustering results, we identify the genes that drive separation between clusters. These marker genes allow us to assign biological meaning to each cluster based on their functional annotation. In the most obvious case, the marker genes for each cluster are a priori associated with particular cell types, allowing us to treat the clustering as a proxy for cell type identity. The same principle can be applied to discover more subtle differences between clusters (e.g., changes in activation or differentiation state) based on the behavior of genes in the affected pathways.

Identification of marker genes is usually based around the retrospective detection of differential expression between clusters. Genes that are more strongly DE are more likely to have caused separate clustering of cells in the first place. Several different statistical tests are available to quantify the differences in expression profiles, and different approaches can be used to consolidate test results into a single ranking of genes for each cluster.

## Learning objectives

* identify genes that differentially expressed between clusters,
* exclusively or not,
* using different methods that test:
  * the mean expression level,
  * the whole distribution,
  * or the proportion of cells expressing the gene
* compile a summary table.

## Load data

We will load the R file keeping the SCE (SingleCellExperiment) object with the outcome of the clustering performed after feature selection and dimensionality reduction of the normalised counts for 500 cells per sample.


```r
# read uncorrected counts
fn <- sprintf("%s/%s/Robjects/%s_sce_nz_postDeconv%s%s_%s_uncorr.Rds",
              projDir,
              outDirBit,
              setName,
              setSuf,
              dsiSuf,
              splSetToGet2) # 'dsi' for data set integration
uncorrected <- readRDS(file=fn)

# Read object in:
#tmpFn <- sprintf("%s/%s/Robjects/%s_sce_nz_postDeconv%s_clustered.Rds", projDir, outDirBit, setName, setSuf)
#tmpFn <- sprintf("%s/%s/Robjects/%s_sce_nz_postDeconv%s_clust.Rds",
#                 projDir, outDirBit, setName, setSuf)
tmpFn <- sprintf("%s/%s/Robjects/%s_sce_nz_postDeconv%s%s_%s_clust.Rds",
                 projDir, outDirBit, setName, setSuf, dsiSuf, splSetToGet2)

print(tmpFn)
```

```
## [1] "/ssd/personal/baller01/20200511_FernandesM_ME_crukBiSs2020/AnaWiSce/AnaKmWiC/Robjects/caron_sce_nz_postDeconv_5hCellPerSpl_dsi_PBMMC_ETV6-RUNX1_clust.Rds"
```

```r
if(!file.exists(tmpFn))
{
	knitr::knit_exit()
}
mnn.out <- readRDS(tmpFn)
rm(tmpFn)

mnn.out <- runUMAP(mnn.out, dimred="corrected")

# copy clustering output to uncorrected SCE
sce <- uncorrected
x <- colData(mnn.out)[,c("Barcode","leiden")] %>% data.frame()
colData(sce) <- colData(uncorrected) %>%
  data.frame() %>%
  dplyr::left_join(x, by="Barcode") %>%
  DataFrame
colData(sce)$cluster <- colData(sce)$leiden
sce$clusterStg <- factor(paste0("c", sce$cluster),
			 levels = paste0("c", levels(sce$cluster)))
rm(uncorrected)
```

Number of cells: .

Number of genes: 16629.

## Detecting genes differentially expressed between clusters

The TSNE plots below show the structure observed with and without data set integration
Cells colored by cluster:


```r
tsney <- plotTSNE(mnn.out, colour_by="leiden") + fontsize
tsnex <- plotTSNE(sce, colour_by="clusterStg") + fontsize
gridExtra::grid.arrange(tsney, tsnex, ncol=2)
```

<img src="clusterMarkerGenes_files/figure-html/unnamed-chunk-2-1.png" width="672" />

Cells colored by sample name:


```r
p1 <- plotTSNE(mnn.out, colour_by="Sample.Name") + fontsize
p2 <- plotTSNE(sce, colour_by="Sample.Name") + fontsize
gridExtra::grid.arrange(p1, p2, ncol=2)
```

<img src="clusterMarkerGenes_files/figure-html/unnamed-chunk-3-1.png" width="672" />

```r
rm(p1, p2)
```

The UMAP plots below show the structure observed with and without data set integration:

Cells colored by cluster:


```r
umapy <- plotUMAP(mnn.out, colour_by="leiden") + fontsize
umapx <- plotUMAP(sce, colour_by="clusterStg") + fontsize
gridExtra::grid.arrange(umapy, umapx, ncol=2)
```

<img src="clusterMarkerGenes_files/figure-html/unnamed-chunk-4-1.png" width="672" />

Cells colored by sample name:


```r
p1 <- plotUMAP(mnn.out, colour_by="Sample.Name") + fontsize
p2 <- plotUMAP(sce, colour_by="Sample.Name") + fontsize
gridExtra::grid.arrange(p1, p2, ncol=2)
```

<img src="clusterMarkerGenes_files/figure-html/unnamed-chunk-5-1.png" width="672" />

```r
rm(p1, p2)
```

### Differential expression analysis

For each cluster, we will identify genes whose expression differ to that of other clusters, for each pair of cluster, using `scran::findMarkers()`. The function  fits a linear model to the log-expression values for each gene using limma [@doi:10.1093/nar/gkv007] and allows testing for differential expression in each cluster compared to the others while accounting for known, uninteresting factors.

We will first identify genes whose average expression differ between clusters, using the Welch t-test (default) with a log-fold change
threshold of 0 (default) and ranking genes based on the outcome of any of the pairwise comparisons (default).


```r
# Clusters called with igraph::cluster_walktrap() are named with digits.
# We add a 'c' prefix (for-cluster) to avoid any confusion: these values are labels, not numerical values

# run scran::findMarkers()
# with default parameters for now
# check the function's manual for details (?scran::findMarkers)
# (the test.type and pval.type options are covered below)
markersWoBatch <- findMarkers(x=sce,
                       groups=sce$clusterStg)
```


```r
markersWiBatch <- findMarkers(x=sce,
                     groups=sce$clusterStg,
                     block=sce$block,
                     #direction="up",
                     #lfc=1,
                     row.data=rowData(sce)[,
                                           #c("Symbol", "ensembl_gene_id"),
                                           ,
                                           drop=FALSE])

 # A cell type of some sort:
demo <- markersWiBatch[["c4"]]
```



```r
as.data.frame(demo[1:20,c("Symbol", "Top", "p.value", "FDR")]) 
```

```
##                 Symbol Top       p.value           FDR
## ENSG00000117632  STMN1   1 1.079294e-262 3.589517e-259
## ENSG00000123416 TUBA1B   1 7.555766e-241 1.794926e-237
## ENSG00000124795    DEK   1  6.490005e-99  1.254910e-96
## ENSG00000164104  HMGB2   1 7.316822e-200 7.604465e-197
## ENSG00000177954  RPS27   1 1.667516e-299 1.386456e-295
## ENSG00000189403  HMGB1   1 5.448322e-223 7.550012e-220
## ENSG00000244734    HBB   1  0.000000e+00  0.000000e+00
## ENSG00000124766   SOX4   2  7.604466e-77  9.802688e-75
## ENSG00000196230   TUBB   2 1.605728e-235 3.337707e-232
## ENSG00000206172   HBA1   2 4.353300e-288 2.413034e-284
## ENSG00000128322  IGLL1   3  2.129167e-89  3.505536e-87
## ENSG00000188536   HBA2   3 1.114332e-285 4.632558e-282
## ENSG00000128218 VPREB3   4 8.331210e-124 2.613956e-121
## ENSG00000164032  H2AFZ   4 5.890852e-158 3.377896e-155
## ENSG00000169877   AHSP   4 1.610936e-241 4.464710e-238
## ENSG00000272398   CD24   4 2.185315e-112 5.118254e-110
## ENSG00000087086    FTL   5 1.168331e-117 3.184947e-115
## ENSG00000122026  RPL21   5 2.077203e-231 3.837980e-228
## ENSG00000166710    B2M   5 1.080283e-207 1.197602e-204
## ENSG00000196549    MME   5  3.392406e-67  3.318371e-65
```

```r
featx <- demo[1,"ensembl_gene_id"]
featy <- demo %>% data.frame %>%
  filter(ensembl_gene_id == featx) %>%
  pull("Symbol")
```


```r
plotExpression(sce,
               x=I(factor(sce$clusterStg)),
               features=featx, # "ENSG00000007312",
               colour_by="block") +
  facet_wrap(~colour_by) +
  ggtitle(glue('{featx} aka {featy}'))
```

<img src="clusterMarkerGenes_files/figure-html/unnamed-chunk-7-1.png" width="672" />

```r
markers <- markersWiBatch
```

Results are compiled in a single table per cluster that stores the outcome of comparisons against the other clusters.
One can then select differentially expressed genes from each pairwise comparison between clusters.

We will define a set of genes for cluster 1 by selecting the top 10 genes of each comparison, and check the test output, eg adjusted p-values and log-fold changes.


```r
clux <- "c1"

# get output table for cluster of interest:
marker.set <- markers[[clux]]
head(marker.set, 3)
```

```
## DataFrame with 3 rows and 23 columns
##                 ensembl_gene_id external_gene_name chromosome_name
##                     <character>        <character>     <character>
## ENSG00000123416 ENSG00000123416             TUBA1B              12
## ENSG00000124766 ENSG00000124766               SOX4               6
## ENSG00000128218 ENSG00000128218             VPREB3              22
##                 start_position end_position strandNum      Symbol
##                      <integer>    <integer> <integer> <character>
## ENSG00000123416       49127782     49131397        -1      TUBA1B
## ENSG00000124766       21593751     21598619         1        SOX4
## ENSG00000128218       23752743     23754425        -1      VPREB3
##                            Type      mean  detected gene_sparsity       Top
##                     <character> <numeric> <numeric>     <numeric> <integer>
## ENSG00000123416 Gene Expression   6.13696   55.0540      0.331165         1
## ENSG00000124766 Gene Expression   2.45547   40.8965      0.387268         1
## ENSG00000128218 Gene Expression   1.11073   28.6456      0.547111         1
##                      p.value          FDR summary.logFC  logFC.c2  logFC.c3
##                    <numeric>    <numeric>     <numeric> <numeric> <numeric>
## ENSG00000123416 1.65043e-162 2.28708e-159      -1.91797  0.719778 0.3009426
## ENSG00000124766  3.84668e-90  1.16302e-87       1.36547  1.077273 0.3147620
## ENSG00000128218 8.91659e-183 2.47123e-179       1.33706  1.337060 0.0150614
##                  logFC.c4  logFC.c5   logFC.c6  logFC.c7  logFC.c8  logFC.c9
##                 <numeric> <numeric>  <numeric> <numeric> <numeric> <numeric>
## ENSG00000123416 -1.917974 -0.468056 -0.0390704  0.577349  0.676016  0.375056
## ENSG00000124766  0.404482  1.192600  0.7622680  1.234901  1.264353  1.365467
## ENSG00000128218  0.174119  1.081052  0.7889714  0.714617  1.302928  1.162344
```

```r
# add gene annotation:
##tmpDf <- marker.set
##tmpDf$ensembl_gene_id <- rownames(tmpDf)
##tmpDf2 <- base::merge(tmpDf, rowData(sce), by="ensembl_gene_id", all.x=TRUE, all.y=F, sort=F)
```

Write table to file:


```r
tmpFn <- sprintf("%s/%s/Robjects/%s_sce_nz_postDeconv%s%s_%s_%s.tsv",
                 projDir, outDirBit, setName, setSuf, dsiSuf, splSetToGet2, clux)
print(tmpFn)
```

```
## [1] "/ssd/personal/baller01/20200511_FernandesM_ME_crukBiSs2020/AnaWiSce/AnaKmWiC/Robjects/caron_sce_nz_postDeconv_5hCellPerSpl_dsi_PBMMC_ETV6-RUNX1_c1.tsv"
```

```r
write.table(marker.set, file=tmpFn, sep="\t", quote=FALSE, row.names=FALSE)
rm(tmpFn)
```

Show expression of marker on t-SNE and UMAP:


```r
mrkSet <- marker.set %>%
  data.frame() %>%
  filter(Symbol %in% rowData(mnn.out)$Symbol) %>%
  arrange(FDR)

geneInd <- 2
geneEns <- mrkSet[geneInd,"ensembl_gene_id"]
geneExt <- mrkSet[geneInd,"external_gene_name"] # external_gene_name ~ Symbol

#tsne1 <- plotTSNE(sce, colour_by=mrkSet[geneInd,"ensembl_gene_id"]) + fontsize
tsne1 <- plotTSNE(sce, colour_by=geneEns) +
  fontsize +
  scale_color_continuous(type = 'viridis',
                         guide = guide_legend(title = geneExt))

umap1 <- plotUMAP(sce, colour_by=geneEns) +
  fontsize +
  scale_color_continuous(type = 'viridis',
                         guide = guide_legend(title = geneExt))

tsne1b <- plotTSNE(mnn.out,
                   colour_by=geneExt,
                   by_exprs_values = "reconstructed") + fontsize

umap1b <- plotUMAP(mnn.out,
                   colour_by=mrkSet[geneInd,"Symbol"],
                   by_exprs_values = "reconstructed") + fontsize
```

Expresion levels of ENSG00000244734 aka HBB:


```r
plotExpression(sce,
               x=I(factor(sce$clusterStg)),
               features=geneEns, # "ENSG00000007312",
               colour_by="block") +
  facet_wrap(~colour_by) +
  ggtitle(glue('{geneEns} aka {geneExt}'))
```

<img src="clusterMarkerGenes_files/figure-html/unnamed-chunk-9-1.png" width="672" />

Without data set integration (sce) 


```r
#tsne1
#umap1
gridExtra::grid.arrange(tsne1, umap1, ncol=2)
```

<img src="clusterMarkerGenes_files/figure-html/unnamed-chunk-10-1.png" width="672" />
With data set integration (mnn.out):


```r
#tsne1b
#umap1b
gridExtra::grid.arrange(tsne1b, umap1b, ncol=2)
```

<img src="clusterMarkerGenes_files/figure-html/unnamed-chunk-11-1.png" width="672" />

With data integration (mnn.out), TSNE with leiden clusters and HBB levels:


```r
gridExtra::grid.arrange(tsney, tsne1b, ncol=2)
```

<img src="clusterMarkerGenes_files/figure-html/unnamed-chunk-12-1.png" width="672" />

Without data integration (sce), TSNE with leiden clusters and HBB levels:


```r
gridExtra::grid.arrange(tsnex, tsne1, ncol=2)
```

<img src="clusterMarkerGenes_files/figure-html/unnamed-chunk-13-1.png" width="672" />

With data integration (mnn.out), UMAP with leiden clusters and HBB levels:


```r
gridExtra::grid.arrange(umapy, umap1b, ncol=2)
```

<img src="clusterMarkerGenes_files/figure-html/unnamed-chunk-14-1.png" width="672" />

Without data integration (sce), UMAP with leiden clusters and HBB levels:


```r
gridExtra::grid.arrange(umapx, umap1, ncol=2)
```

<img src="clusterMarkerGenes_files/figure-html/unnamed-chunk-15-1.png" width="672" />


```r
rm(tsne1, tsne1b, umap1, umap1b)
rm(tsnex, tsney, umapx, umapy)
```

Gene set enrichment analyses used for bulk RNA-seq may be used to characterise clusters further. 

### Heatmap

As for bulk RNA, differences in expression profiles of the top genes can be visualised with a heatmap. 

Normalised counts:


```r
# select some top genes:
top.markers <- rownames(marker.set)[marker.set$Top <= 10]

# have matrix to annotate sample with cluster and sample:
tmpData <- logcounts(sce)[top.markers,]
# concat sample and barcode names to make unique name across the whole data set
tmpCellNames <- paste(colData(sce)$Sample.Name, colData(sce)$Barcode, sep="_")
# use these to namecolumn of matrix the show as heatmap:
colnames(tmpData) <- tmpCellNames # colData(sce)$Barcode                    

# columns annotation with cell name:
mat_col <- data.frame(cluster = sce$cluster,
		      sample = sce$Sample.Name,
		      type = sce$source_name
		)
rownames(mat_col) <- colnames(tmpData)
rownames(mat_col) <- tmpCellNames # colData(sce)$Barcode

# Prepare colours for clusters:
colourCount = length(unique(sce$cluster))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

mat_colors <- list(group = getPalette(colourCount))
names(mat_colors$group) <- unique(sce$cluster)

breaksVec = seq(-5, 5, by = 0.1)

# plot heatmap:
pheatmap(tmpData,
           border_color      = NA,
           show_colnames     = FALSE,
           #show_rownames     = FALSE,
           show_rownames     = TRUE,
           drop_levels       = TRUE,
           labels_row        = rowData(sce)[rownames(tmpData),"Symbol"],
           annotation_col    = mat_col,
           annotation_colors = mat_colors,
           color             = colorRampPalette(
             rev(brewer.pal(n = 7,
                            name = "RdYlBu")))(length(breaksVec)),
           breaks            = breaksVec,
	   fontsize_row      = 7
           )
```

<img src="clusterMarkerGenes_files/figure-html/marker_set_clu1_heatmap_unsorted-1.png" width="672" />


```r
# remove MALAT1
tmpData2 <- tmpData %>%
  as.matrix() %>%
  data.frame() %>%
  tibble::rownames_to_column("ensId") %>%
  filter(ensId != "ENSG00000251562") %>%
  tibble::column_to_rownames("ensId")
colnames(tmpData2) <- gsub("\\.", "-", colnames(tmpData2))

pheatmap(tmpData2,
           border_color      = NA,
           show_colnames     = FALSE,
           #show_rownames     = FALSE,
           show_rownames     = TRUE,
           drop_levels       = TRUE,
           labels_row        = rowData(sce)[rownames(tmpData2),"Symbol"],
           annotation_col    = mat_col,
           annotation_colors = mat_colors,
           color             = colorRampPalette(
             rev(brewer.pal(n = 7,
                            name = "RdYlBu")))(length(breaksVec)),
           breaks            = breaksVec,
	   fontsize_row      = 7
           )
rm(tmpData2)
```

<!--
One can sort both the gene and sample dendrograms to improve the heatmap.
-->


```r
library(dendsort)

mat <- tmpData
mat_cluster_cols <- hclust(dist(t(mat)))

sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

mat_cluster_cols <- sort_hclust(mat_cluster_cols)
#plot(mat_cluster_cols, main = "Sorted Dendrogram", xlab = "", sub = "")

mat_cluster_rows <- sort_hclust(hclust(dist(mat)))
rm(mat)

pheatmap(tmpData,
           border_color      = NA,
           show_colnames     = FALSE,
           show_rownames     = FALSE,
           drop_levels       = TRUE,
           annotation_col    = mat_col,
           annotation_colors = mat_colors,
           cluster_cols      = mat_cluster_cols,
           cluster_rows      = mat_cluster_rows
         )
```

Counts corrected for batch effect ('corrected'):


```r
# select some top genes:
#top.markers <- rownames(marker.set)[marker.set$Top <= 10]
top.markers <- marker.set[marker.set$Top <= 10,] %>%
  data.frame() %>%
  pull(Symbol)

# have matrix to annotate sample with cluster and sample:
tmpData3 <- assay(mnn.out) %>%
  as.matrix() %>%
  data.frame() %>%
  tibble::rownames_to_column("Symbol") %>%
  filter(Symbol %in% top.markers) %>%
  tibble::column_to_rownames("Symbol")
colnames(tmpData3) <- colData(mnn.out)$Barcode
colnames(tmpData3) <- gsub("\\.", "-", colnames(tmpData3))

# concat sample and barcode names to make unique name across the whole data set
tmpCellNames <- paste(colData(mnn.out)$Sample.Name,
                      colData(mnn.out)$Barcode, sep="_")
# use these to name columns of matrix the show as heatmap:
colnames(tmpData3) <- tmpCellNames # colData(sce)$Barcode                    

# columns annotation with cell name:
mat_col <- data.frame(cluster = sce$cluster,
		      sample = sce$Sample.Name,
		      type = sce$source_name
		)
rownames(mat_col) <- colnames(tmpData)
rownames(mat_col) <- tmpCellNames # colData(sce)$Barcode

# Prepare colours for clusters:
colourCount = length(unique(sce$cluster))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

mat_colors <- list(group = getPalette(colourCount))
names(mat_colors$group) <- unique(sce$cluster)

# plot heatmap:
pheatmap(tmpData3[,order(colData(sce)$clusterStg)],
           border_color      = NA,
           show_colnames     = FALSE,
           #show_rownames     = FALSE,
           show_rownames     = TRUE,
           drop_levels       = TRUE,
           labels_row        = rowData(sce)[rownames(tmpData),"Symbol"],
           annotation_col    = mat_col,
           annotation_colors = mat_colors,
           cluster_cols      = FALSE,
           fontsize_row      = 7
           )
```

<img src="clusterMarkerGenes_files/figure-html/marker_set_clu1_heatmap_unsorted_recon-1.png" width="672" />

```r
rm(tmpData3)
rm(mnn.out, tmpData)
```

To demonstrate how to interpret the results, we will use cluster 9 as our cluster of interest. The relevant DataFrame contains log2-fold changes of expression in cluster 9 over each other cluster, along with several statistics obtained by combining p-values (Simes 1986) across the pairwise comparisons involving cluster 9.


```r
chosen <- "c9"
interesting <- markers[[chosen]]
print(colnames(interesting))
```

```
##  [1] "ensembl_gene_id"    "external_gene_name" "chromosome_name"   
##  [4] "start_position"     "end_position"       "strandNum"         
##  [7] "Symbol"             "Type"               "mean"              
## [10] "detected"           "gene_sparsity"      "Top"               
## [13] "p.value"            "FDR"                "summary.logFC"     
## [16] "logFC.c1"           "logFC.c2"           "logFC.c3"          
## [19] "logFC.c4"           "logFC.c5"           "logFC.c6"          
## [22] "logFC.c7"           "logFC.c8"
```

Of particular interest is the Top field. The set of genes with Top ≤X is the union of the top X genes (ranked by p-value) from each pairwise comparison involving cluster 9. For example, the set of all genes with Top values of 1 contains the gene with the lowest p-value from each comparison. Similarly, the set of genes with Top values less than or equal to 10 contains the top 10 genes from each comparison. The Top field represents findMarkers()’s approach to consolidating multiple pairwise comparisons into a single ranking for each cluster; each DataFrame produced by findMarkers() will order genes based on the Top value by default.


```r
interesting.coluToShow <- c("Symbol", "Top", "FDR", "summary.logFC")
interesting[1:10,interesting.coluToShow]
```

```
## DataFrame with 10 rows and 4 columns
##                      Symbol       Top          FDR summary.logFC
##                 <character> <integer>    <numeric>     <numeric>
## ENSG00000090013       BLVRB         1 5.73724e-151     -4.045409
## ENSG00000124766        SOX4         1  7.99580e-87     -1.365467
## ENSG00000124795         DEK         1  1.79870e-95     -1.053542
## ENSG00000156738       MS4A1         1  4.06226e-41     -1.557764
## ENSG00000204482        LST1         1  3.17916e-52     -1.154475
## ENSG00000227507         LTB         1  1.52234e-81     -1.419390
## ENSG00000271503        CCL5         1  5.71426e-27     -2.236276
## ENSG00000004939      SLC4A1         2 1.49020e-122     -3.175878
## ENSG00000082074        FYB1         2  4.57767e-74     -0.882386
## ENSG00000105374        NKG7         2  2.62740e-26     -2.054456
```

We use the Top field to identify a set of genes that is guaranteed to distinguish cluster 9 from any other cluster. Here, we examine the top 6 genes from each pairwise comparison.


```r
best.set <- interesting[interesting$Top <= 6,]
logFCs <- getMarkerEffects(best.set)
logFCs.ens <- rownames(logFCs)
rownames(logFCs) <- rowData(sce)[rownames(logFCs), "Symbol"]

library(pheatmap)
pheatmap(logFCs, breaks=seq(-5, 5, length.out=101))
```

<img src="clusterMarkerGenes_files/figure-html/unnamed-chunk-19-1.png" width="576" />

## Using the log-fold change

Our previous findMarkers() call considers both up- and downregulated genes to be potential markers. However, downregulated genes are less appealing as markers as it is more difficult to interpret and experimentally validate an absence of expression. To focus on up-regulated markers, we can instead perform a one-sided t-test to identify genes that are upregulated in each cluster compared to the others. This is achieved by setting direction="up" in the findMarkers() call.


```r
markersWoBatch.up <- findMarkers(sce,
				 groups = sce$clusterStg,
				 direction = "up")
```


```r
markersWiBatch.up <- findMarkers(x = sce,
                     groups = sce$clusterStg,
                     block = sce$block,
                     direction = "up",
                     #lfc=1,
                     row.data = rowData(sce)[,
                                           #c("Symbol", "ensembl_gene_id"),
                                           ,
                                           drop=FALSE])
markers.up <- markersWiBatch.up
```


```r
interesting.up <- markers.up[[chosen]]
interesting.up[1:10,interesting.coluToShow]
```

```
## DataFrame with 10 rows and 4 columns
##                      Symbol       Top         FDR summary.logFC
##                 <character> <integer>   <numeric>     <numeric>
## ENSG00000051523        CYBA         1 1.46243e-09       2.10059
## ENSG00000113387        SUB1         1 2.70422e-09       2.12860
## ENSG00000137818       RPLP1         1 2.70422e-09       1.62362
## ENSG00000166710         B2M         1 2.68397e-08       1.84380
## ENSG00000140988        RPS2         2 6.10877e-09       1.36627
## ENSG00000166562      SEC11C         2 7.40121e-08       2.01498
## ENSG00000170476        MZB1         3 2.69120e-07       3.31396
## ENSG00000180879        SSR4         3 1.75120e-07       2.50935
## ENSG00000234745       HLA-B         3 4.35316e-08       1.89577
## ENSG00000172183       ISG20         4 7.19818e-07       1.35939
```

The t-test also allows us to specify a non-zero log-fold change as the null hypothesis. This allows us to consider the magnitude of the log-fold change in our p-value calculations, in a manner that is more rigorous than simply filtering directly on the log-fold changes (McCarthy and Smyth 2009). (Specifically, a simple threshold does not consider the variance and can enrich for genes that have both large log-fold changes and large variances.) We perform this by setting lfc= in our findMarkers() call - when combined with direction=, this tests for genes with log-fold changes that are significantly greater than 1:


```r
markersWoBatch.up2 <- findMarkers(sce, groups=sce$clusterStg, direction="up", lfc=1)
```


```r
markersWiBatch.up2 <- findMarkers(x=sce,
                     groups=sce$clusterStg,
                     block=sce$block,
                     direction="up",
                     #lfc=1,
                     row.data=rowData(sce)[,
                                           #c("Symbol", "ensembl_gene_id"),
                                           ,
                                           drop=FALSE])
markers.up2 <- markersWiBatch.up2
```


```r
interesting.up2 <- markers.up2[[chosen]]
interesting.up2[1:10,interesting.coluToShow]
```

```
## DataFrame with 10 rows and 4 columns
##                      Symbol       Top         FDR summary.logFC
##                 <character> <integer>   <numeric>     <numeric>
## ENSG00000051523        CYBA         1 1.46243e-09       2.10059
## ENSG00000113387        SUB1         1 2.70422e-09       2.12860
## ENSG00000137818       RPLP1         1 2.70422e-09       1.62362
## ENSG00000166710         B2M         1 2.68397e-08       1.84380
## ENSG00000140988        RPS2         2 6.10877e-09       1.36627
## ENSG00000166562      SEC11C         2 7.40121e-08       2.01498
## ENSG00000170476        MZB1         3 2.69120e-07       3.31396
## ENSG00000180879        SSR4         3 1.75120e-07       2.50935
## ENSG00000234745       HLA-B         3 4.35316e-08       1.89577
## ENSG00000172183       ISG20         4 7.19818e-07       1.35939
```

These two settings yield a more focused set of candidate marker genes that are upregulated in cluster 9.


```r
best.set <- interesting.up2[interesting.up2$Top <= 5,]
logFCs <- getMarkerEffects(best.set)
logFCs.ens <- rownames(logFCs)
rownames(logFCs) <- rowData(sce)[rownames(logFCs), "Symbol"]

library(pheatmap)
pheatmap(logFCs, breaks=seq(-5, 5, length.out=101))
```

<img src="clusterMarkerGenes_files/figure-html/unnamed-chunk-26-1.png" width="672" />

Of course, this increased stringency is not without cost. If only upregulated genes are requested from findMarkers(), any cluster defined by downregulation of a marker gene will not contain that gene among the top set of features in its DataFrame. This is occasionally relevant for subtypes or other states that are distinguished by high versus low expression of particular genes. Similarly, setting an excessively high log-fold change threshold may discard otherwise useful genes. For example, a gene upregulated in a small proportion of cells of a cluster will have a small log-fold change but can still be an effective marker if the focus is on specificity rather than sensitivity.

## Finding cluster-specific markers

By default, findMarkers() will give a high ranking to genes that are differentially expressed in any pairwise comparison. This is because a gene only needs a very low p
-value in a single pairwise comparison to achieve a low Top value. A more stringent approach would only consider genes that are differentially expressed in all pairwise comparisons involving the cluster of interest. To achieve this, we set pval.type="all" in findMarkers() to use an intersection-union test (Berger and Hsu 1996) where the combined p-value for each gene is the maximum of the p-values from all pairwise comparisons. A gene will only achieve a low combined p-value if it is strongly DE in all comparisons to other clusters.


```r
# We can combine this with 'direction='.
markersWoBatch.up3 <- findMarkers(sce,
				  groups=sce$clusterStg,
				  pval.type="all",
				  direction="up")
```


```r
markersWiBatch.up3 <- findMarkers(x=sce,
                     groups=sce$clusterStg,
                     block=sce$block,
                     pval.type="all",
                     direction="up",
                     #lfc=1,
                     row.data=rowData(sce)[,
                                           #c("Symbol", "ensembl_gene_id"),
                                           ,
                                           drop=FALSE])
markers.up3 <- markersWiBatch.up3
```


```r
interesting.up3 <- markers.up3[[chosen]]
interesting.colu <- intersect(interesting.coluToShow, colnames(interesting.up3))
interesting.up3[1:10,interesting.colu]
```

```
## DataFrame with 10 rows and 3 columns
##                      Symbol         FDR summary.logFC
##                 <character>   <numeric>     <numeric>
## ENSG00000113387        SUB1 2.59042e-07       1.88802
## ENSG00000166562      SEC11C 2.59042e-07       1.91824
## ENSG00000180879        SSR4 2.38922e-06       2.37547
## ENSG00000132465      JCHAIN 4.57138e-06       4.96266
## ENSG00000134285      FKBP11 2.35844e-05       1.71931
## ENSG00000048462    TNFRSF17 2.35844e-05       1.60136
## ENSG00000166794        PPIB 2.35844e-05       1.87611
## ENSG00000170476        MZB1 2.35844e-05       2.31249
## ENSG00000166598     HSP90B1 2.72866e-05       1.96988
## ENSG00000106803      SEC61B 1.64771e-04       1.20282
```

This strategy will only report genes that are highly specific to the cluster of interest. When it works, it can be highly effective as it generates a small focused set of candidate markers. However, any gene that is expressed at the same level in two or more clusters will simply not be detected. This is likely to discard many interesting genes, especially if the clusters are finely resolved with weak separation. To give a concrete example, consider a mixed population of CD4+-only, CD8+-only, double-positive and double-negative T cells. With pval.type="all", neither Cd4 or Cd8 would be detected as subpopulation-specific markers because each gene is expressed in two subpopulations. In comparison, pval.type="any" will detect both of these genes as they will be DE between at least one pair of subpopulations.

If pval.type="all" is too stringent yet pval.type="any" is too generous, a compromise is to set pval.type="some". For each gene, we apply the Holm-Bonferroni correction across its p
-values and take the middle-most value as the combined p-value. This effectively tests the global null hypothesis that at least 50% of the individual pairwise comparisons exhibit no DE. We then rank the genes by their combined p-values to obtain an ordered set of marker candidates. The aim is to improve the conciseness of the top markers for defining a cluster while mitigating the risk of discarding useful genes that are not DE to all other clusters. The downside is that taking this compromise position sacrifices the theoretical guarantees offered at the other two extremes.


```r
markersWoBatch.up4 <- findMarkers(sce,
				  groups=sce$clusterStg,
				  pval.type="some",
				  direction="up")
```


```r
markersWiBatch.up4 <- findMarkers(x=sce,
                     groups=sce$clusterStg,
                     block=sce$block,
                     pval.type="some",
                     direction="up",
                     #lfc=1,
                     row.data=rowData(sce)[,
                                           #c("Symbol", "ensembl_gene_id"),
                                           ,
                                           drop=FALSE])
```


```r
markers.up4 <- markersWiBatch.up4
interesting.up4 <- markers.up4[[chosen]]
interesting.colu <- intersect(interesting.coluToShow, colnames(interesting.up4))
interesting.up4[1:10,interesting.colu]
```

```
## DataFrame with 10 rows and 3 columns
##                      Symbol         FDR summary.logFC
##                 <character>   <numeric>     <numeric>
## ENSG00000113387        SUB1 5.95797e-08       2.02728
## ENSG00000051523        CYBA 3.93001e-07       1.76665
## ENSG00000166562      SEC11C 6.15867e-07       1.95627
## ENSG00000180879        SSR4 9.31346e-07       2.60317
## ENSG00000170476        MZB1 1.43640e-06       2.74844
## ENSG00000132465      JCHAIN 5.53928e-06       5.41802
## ENSG00000166598     HSP90B1 2.03107e-05       2.13393
## ENSG00000166794        PPIB 2.03107e-05       1.87550
## ENSG00000172183       ISG20 2.03107e-05       1.12174
## ENSG00000134285      FKBP11 2.44952e-05       1.91535
```

In both cases, a different method is used to compute the summary effect size compared to pval.type="any". For pval.type="all", the summary log-fold change is defined as that corresponding to the pairwise comparison with the largest p-value, while for pval.type="some", it is defined as the log-fold change for the comparison with the middle-most p-value. This reflects the calculation of the combined p-value and avoids focusing on genes with strong changes in only one comparison.

## Using the Wilcoxon rank sum test

The Wilcoxon rank sum test (also known as the Wilcoxon-Mann-Whitney test, or WMW test) is another widely used method for pairwise comparisons between groups of observations. Its strength lies in the fact that it directly assesses separation between the expression distributions of different clusters. The WMW test statistic is proportional to the area-under-the-curve (AUC), i.e., the concordance probability, which is the probability of a random cell from one cluster having higher expression than a random cell from another cluster. In a pairwise comparison, AUCs of 1 or 0 indicate that the two clusters have perfectly separated expression distributions. Thus, the WMW test directly addresses the most desirable property of a candidate marker gene, while the t-test only does so indirectly via the difference in the means and the intra-group variance.

We perform WMW tests by again using the findMarkers() function, this time with test="wilcox". This returns a list of DataFrames containing ranked candidate markers for each cluster. The direction=, lfc= and pval.type= arguments can be specified and have the same interpretation as described for t-tests. We demonstrate below by detecting upregulated genes in each cluster with direction="up".


```r
markers.wmw <- findMarkers(sce, groups=sce$clusterStg, test="wilcox", direction="up")
```


```r
markersWiBatch.wmw <- findMarkers(x=sce,
                     groups=sce$clusterStg,
                     block=sce$block,
                     test="wilcox",
                     direction="up",
                     row.data=rowData(sce)[,
                                           #c("Symbol", "ensembl_gene_id"),
                                           ,
                                           drop=FALSE])
```


```r
markers.wmw <- markersWiBatch.wmw
#print(names(markers.wmw))
```

To explore the results in more detail, we focus on the DataFrame for cluster 9. The interpretation of Top is the same as described for t-tests, and Simes’ method is again used to combine p-values across pairwise comparisons. If we want more focused sets, we can also change pval.type= as previously described.


```r
interesting.coluToShow <- c("Symbol", "Top", "FDR", "summary.AUC")

interesting.wmw <- markers.wmw[[chosen]]
interesting.wmw[1:10,interesting.coluToShow]
```

```
## DataFrame with 10 rows and 4 columns
##                      Symbol       Top          FDR summary.AUC
##                 <character> <integer>    <numeric>   <numeric>
## ENSG00000048462    TNFRSF17         1 3.68222e-123    0.996267
## ENSG00000134285      FKBP11         1  9.56095e-83    0.966667
## ENSG00000254709       IGLL5         1 5.05036e-161    0.979676
## ENSG00000070081       NUCB2         2  1.49469e-26    0.925664
## ENSG00000099958       DERL3         2  1.30411e-84    0.917676
## ENSG00000110777     POU2AF1         2  3.24097e-68    0.888594
## ENSG00000138185      ENTPD1         2 1.22295e-133    0.895064
## ENSG00000165272        AQP3         2 3.14631e-107    0.972418
## ENSG00000103811        CTSH         3  2.08541e-73    0.851989
## ENSG00000124772       CPNE5         3  4.73626e-63    0.788544
```

The DataFrame contains the AUCs from comparing cluster 9 to every other cluster. A value greater than 0.5 indicates that the gene is upregulated in the current cluster compared to the other cluster, while values less than 0.5 correspond to downregulation. We would typically expect AUCs of 0.7-0.8 for a strongly upregulated candidate marker.


```r
best.set <- interesting.wmw[interesting.wmw$Top <= 5,]
AUCs <- getMarkerEffects(best.set, prefix="AUC")
AUCs.ens <- rownames(AUCs)
rownames(AUCs) <- rowData(sce)[rownames(AUCs), "Symbol"]


library(pheatmap)
pheatmap(AUCs,
	 breaks = seq(0, 1, length.out=21),
	 color = viridis::viridis(21))
```

<img src="clusterMarkerGenes_files/figure-html/unnamed-chunk-37-1.png" width="672" />

One practical advantage of the WMW test over the Welch t-test is that it is symmetric with respect to differences in the size of the groups being compared. This means that, all else being equal, the top-ranked genes on each side of a DE comparison will have similar expression profiles regardless of the number of cells in each group. In contrast, the t-test will favor genes where the larger group has the higher relative variance as this increases the estimated degrees of freedom and decreases the resulting p-value. This can lead to unappealing rankings when the aim is to identify genes upregulated in smaller groups. The WMW test is not completely immune to variance effects - for example, it will slightly favor detection of DEGs at low average abundance where the greater number of ties at zero deflates the approximate variance of the rank sum statistic - but this is relatively benign as the selected genes are still fairly interesting.

<!-- We observe both of these effects in a comparison between alpha and gamma cells in the human pancreas data set from Lawlor et al. (2017) (Figure 11.4). -->


```r
print(head(sce$clusterStg))
```


```r
type1 <- "c1"
#type2 <- "c10"
type2 <- levels(sce$clusterStg)[nlevels(sce$clusterStg)]
marker.t <- findMarkers(sce,
                        groups=sce$clusterStg,
                        direction="up",
                        block=sce$block,
                        restrict=c(type1, type2),
                        row.data=rowData(sce)[,
                                           c("Symbol", "ensembl_gene_id"),
                                           drop=FALSE])
marker.w <- findMarkers(sce,
                        groups=sce$clusterStg,
                        direction="up",
                        block=sce$block,
                        restrict=c(type1, type2),
                        test.type="wilcox",
                        row.data=rowData(sce)[,
                                           c("Symbol", "ensembl_gene_id"),
                                           drop=FALSE]
                        )
```


```r
# Upregulated in type 1:
marker.type1.t <- marker.t[[type1]]
marker.type1.w <- marker.w[[type1]]
chosen.type1.t <- rownames(marker.type1.t)[1:30]
chosen.type1.w <- rownames(marker.type1.w)[1:30]
u.type1.t <- setdiff(chosen.type1.t, chosen.type1.w)
u.type1.w <- setdiff(chosen.type1.w, chosen.type1.t)

# Upregulated in gamma:
marker.type2.t <- marker.t[[type2]]
marker.type2.w <- marker.w[[type2]]
chosen.type2.t <- rownames(marker.type2.t)[1:30]
chosen.type2.w <- rownames(marker.type2.w)[1:30]
u.type2.t <- setdiff(chosen.type2.t, chosen.type2.w)
u.type2.w <- setdiff(chosen.type2.w, chosen.type2.t)

u.type1.t <- u.type1.t[1:min(4,length(u.type1.t))]
u.type1.w <- u.type1.w[1:min(4,length(u.type1.w))]
u.type2.t <- u.type2.t[1:min(4,length(u.type2.t))]
u.type2.w <- u.type2.w[1:min(4,length(u.type2.w))]
```

u.type1.t


```r
marker.type1.t[,] %>%
  data.frame() %>%
  filter(ensembl_gene_id %in% u.type1.t)
```

```
##                  Symbol ensembl_gene_id Top      p.value          FDR
## ENSG00000242574 HLA-DMB ENSG00000242574   4 7.999069e-43 3.325413e-39
## ENSG00000066923   STAG3 ENSG00000066923   5 2.123082e-40 7.060948e-37
## ENSG00000164330    EBF1 ENSG00000164330   6 1.970585e-39 5.461476e-36
## ENSG00000124795     DEK ENSG00000124795   7 1.384151e-38 3.288151e-35
##                 summary.logFC  logFC.c9
## ENSG00000242574     0.3230307 0.3230307
## ENSG00000066923     0.3309251 0.3309251
## ENSG00000164330     0.4874803 0.4874803
## ENSG00000124795     0.4620890 0.4620890
```

u.type1.w


```r
marker.type1.w[,] %>%
  data.frame() %>%
  filter(ensembl_gene_id %in% u.type1.w)
```

```
##                   Symbol ensembl_gene_id Top      p.value        FDR
## ENSG00000204287  HLA-DRA ENSG00000204287   1 1.981116e-06 0.03294397
## ENSG00000117632    STMN1 ENSG00000117632   2 2.729759e-05 0.09106999
## ENSG00000223865 HLA-DPB1 ENSG00000223865   3 3.244727e-05 0.09106999
## ENSG00000019582     CD74 ENSG00000019582   4 3.474180e-05 0.09106999
##                 summary.AUC    AUC.c9
## ENSG00000204287   0.9933637 0.9933637
## ENSG00000117632   0.9321858 0.9321858
## ENSG00000223865   0.9299046 0.9299046
## ENSG00000019582   0.9311489 0.9311489
```

u.type2.t


```r
marker.type2.t[,] %>%
  data.frame() %>%
  filter(ensembl_gene_id %in% u.type2.t)
```

```
##                 Symbol ensembl_gene_id Top      p.value          FDR
## ENSG00000051523   CYBA ENSG00000051523   1 9.320325e-13 1.257556e-08
## ENSG00000113387   SUB1 ENSG00000113387   2 1.512485e-12 1.257556e-08
## ENSG00000166562 SEC11C ENSG00000166562   3 6.632282e-12 3.676274e-08
## ENSG00000180879   SSR4 ENSG00000180879   4 1.053099e-11 4.377998e-08
##                 summary.logFC logFC.c1
## ENSG00000051523      1.680844 1.680844
## ENSG00000113387      2.107369 2.107369
## ENSG00000166562      2.014979 2.014979
## ENSG00000180879      2.509347 2.509347
```

u.type2.w


```r
marker.type2.w[,] %>%
  data.frame() %>%
  filter(ensembl_gene_id %in% u.type2.w)
```

```
##                    Symbol ensembl_gene_id Top       p.value           FDR
## ENSG00000254709     IGLL5 ENSG00000254709   1 3.796350e-166 6.312951e-162
## ENSG00000138185    ENTPD1 ENSG00000138185   2 2.757867e-138 1.528685e-134
## ENSG00000240505 TNFRSF13B ENSG00000240505   3 2.757867e-138 1.528685e-134
## ENSG00000211890     IGHA2 ENSG00000211890   5 8.148131e-116 2.709905e-112
##                 summary.AUC    AUC.c1
## ENSG00000254709   0.9796765 0.9796765
## ENSG00000138185   0.8950643 0.8950643
## ENSG00000240505   0.8950643 0.8950643
## ENSG00000211890   0.8947532 0.8947532
```


```r
# Examining all uniquely detected markers in each direction.
library(scater)
subset <- sce[,sce$clusterStg %in% c(type1, type2)]
gridExtra::grid.arrange(
    plotExpression(subset, x="clusterStg", features=u.type1.t, ncol=2) +
        ggtitle(sprintf("Upregulated in %s, t-test-only", type1)),
    plotExpression(subset, x="clusterStg", features=u.type1.w, ncol=2) +
        ggtitle(sprintf("Upregulated in %s, WMW-test-only", type1)),
    plotExpression(subset, x="clusterStg", features=u.type2.t, ncol=2) +
        ggtitle(sprintf("Upregulated in %s, t-test-only", type2)),
    plotExpression(subset, x="clusterStg", features=u.type2.w, ncol=2) +
        ggtitle(sprintf("Upregulated in %s, WMW-test-only", type2)),
    ncol=2
)
```

<img src="clusterMarkerGenes_files/figure-html/unnamed-chunk-42-1.png" width="672" />

```r
rm(u.type1.t, u.type1.w, u.type2.t, u.type2.w)
```

The main disadvantage of the WMW test is that the AUCs are much slower to compute compared to t-statistics. This may be inconvenient for interactive analyses involving multiple iterations of marker detection. We can mitigate this to some extent by parallelizing these calculations using the BPPARAM= argument in findMarkers().

##  Using a binomial test

The binomial test identifies genes that differ in the proportion of expressing cells between clusters. (For the purposes of this section, a cell is considered to express a gene simply if it has non-zero expression for that gene.) This represents a much more stringent definition of marker genes compared to the other methods, as differences in expression between clusters are effectively ignored if both distributions of expression values are not near zero. The premise is that genes are more likely to contribute to important biological decisions if they were active in one cluster and silent in another, compared to more subtle “tuning” effects from changing the expression of an active gene. From a practical perspective, a binary measure of presence/absence is easier to validate.

We perform pairwise binomial tests between clusters using the findMarkers() function with test="binom". This returns a list of DataFrames containing marker statistics for each cluster such as the Top rank and its p-value. Here, the effect size is reported as the log-fold change in this proportion between each pair of clusters. Large positive log-fold changes indicate that the gene is more frequently expressed in one cluster compared to the other. We focus on genes that are upregulated in each cluster compared to the others by setting direction="up".


```r
markersWoBatch.binom <- findMarkers(sce,
				    test="binom",
				    direction="up",
				    groups=sce$clusterStg)
```



```r
markersWiBatch.binom <- findMarkers(x=sce,
                     groups=sce$clusterStg,
                     block=sce$block,
                     test="binom",
                     direction="up",
                     row.data=rowData(sce)[,
                                           #c("Symbol", "ensembl_gene_id"),
                                           ,
                                           drop=FALSE])
```


```r
markers.binom <- markersWiBatch.binom
print(names(markers.binom))
```

```
## [1] "c1" "c2" "c3" "c4" "c5" "c6" "c7" "c8" "c9"
```


```r
interesting.binom <- markers.binom[[chosen]]
print(colnames(interesting.binom))
```

```
##  [1] "ensembl_gene_id"    "external_gene_name" "chromosome_name"   
##  [4] "start_position"     "end_position"       "strandNum"         
##  [7] "Symbol"             "Type"               "mean"              
## [10] "detected"           "gene_sparsity"      "Top"               
## [13] "p.value"            "FDR"                "summary.logFC"     
## [16] "logFC.c1"           "logFC.c2"           "logFC.c3"          
## [19] "logFC.c4"           "logFC.c5"           "logFC.c6"          
## [22] "logFC.c7"           "logFC.c8"
```

The plot below confirms that the top genes exhibit strong differences in the proportion of expressing cells in cluster 9 compared to the others.


```r
top.genes <- head(rownames(interesting.binom))
#plotExpression(sce, x="clusterStg", features=top.genes)
plotExpression(sce, x="clusterStg",
               colour_by="clusterStg",
               features=top.genes[1:4] )
```

<img src="clusterMarkerGenes_files/figure-html/unnamed-chunk-47-1.png" width="672" />




```r
plotExpression(sce, x="clusterStg",
               colour_by="clusterStg",
               features=top.genes[1] )
```

<img src="clusterMarkerGenes_files/figure-html/unnamed-chunk-48-1.png" width="672" />

```r
plotExpression(sce, x="clusterStg",
               colour_by="clusterStg",
               features=top.genes[2] )
```

<img src="clusterMarkerGenes_files/figure-html/unnamed-chunk-48-2.png" width="672" />

```r
plotExpression(sce, x="clusterStg",
               colour_by="clusterStg",
               features=top.genes[3] )
```

<img src="clusterMarkerGenes_files/figure-html/unnamed-chunk-48-3.png" width="672" />

```r
plotExpression(sce, x="clusterStg",
               colour_by="clusterStg",
               features=top.genes[4] )
```

<img src="clusterMarkerGenes_files/figure-html/unnamed-chunk-48-4.png" width="672" />

## Combining multiple marker statistics

On occasion, we might want to combine marker statistics from several testing regimes into a single DataFrame. This allows us to easily inspect multiple statistics at once to verify that a particular gene is a strong candidate marker. For example, a large AUC from the WMW test indicates that the expression distributions are well-separated between clusters, while the log-fold change reported with the t-test provides a more interpretable measure of the magnitude of the change in expression. We use the multiMarkerStats() to merge the results of separate findMarkers() calls into one DataFrame per cluster, with statistics interleaved to facilitate a direct comparison between different test regimes.


```r
combined <- multiMarkerStats(
    t=findMarkers(sce, groups=sce$clusterStg, direction="up"),
    wilcox=findMarkers(sce, groups=sce$clusterStg, test="wilcox", direction="up"),
    binom=findMarkers(sce, groups=sce$clusterStg, test="binom", direction="up")
)
```


```r
combined <- multiMarkerStats(
    t=findMarkers(sce, groups=sce$clusterStg, direction="up", block=sce$block),
    wilcox=findMarkers(sce, groups=sce$clusterStg, test="wilcox", direction="up", block=sce$block),
    binom=findMarkers(sce, groups=sce$clusterStg, test="binom", direction="up", block=sce$block)
)
```


```r
# Interleaved marker statistics from both tests for each cluster.
print(colnames(combined[["c1"]]))
```

```
##  [1] "Top"                 "p.value"             "FDR"                
##  [4] "t.Top"               "wilcox.Top"          "binom.Top"          
##  [7] "t.p.value"           "wilcox.p.value"      "binom.p.value"      
## [10] "t.FDR"               "wilcox.FDR"          "binom.FDR"          
## [13] "t.summary.logFC"     "wilcox.summary.AUC"  "binom.summary.logFC"
## [16] "t.logFC.c2"          "wilcox.AUC.c2"       "binom.logFC.c2"     
## [19] "t.logFC.c3"          "wilcox.AUC.c3"       "binom.logFC.c3"     
## [22] "t.logFC.c4"          "wilcox.AUC.c4"       "binom.logFC.c4"     
## [25] "t.logFC.c5"          "wilcox.AUC.c5"       "binom.logFC.c5"     
## [28] "t.logFC.c6"          "wilcox.AUC.c6"       "binom.logFC.c6"     
## [31] "t.logFC.c7"          "wilcox.AUC.c7"       "binom.logFC.c7"     
## [34] "t.logFC.c8"          "wilcox.AUC.c8"       "binom.logFC.c8"     
## [37] "t.logFC.c9"          "wilcox.AUC.c9"       "binom.logFC.c9"
```

```r
#head(combined[["c1"]][,1:9])
combined[["c1"]]$Symbol <- rowData(sce)[rownames(combined[["c1"]]), "Symbol"]
tmpCol <- c("Symbol", colnames(combined[["c1"]])[1:9])
head(combined[["c1"]][,tmpCol])
```

```
## DataFrame with 6 rows and 10 columns
##                      Symbol       Top     p.value         FDR     t.Top
##                 <character> <integer>   <numeric>   <numeric> <integer>
## ENSG00000124766        SOX4         3 1.49550e-41 1.18422e-38         1
## ENSG00000128218      VPREB3         4 6.04516e-64 1.00525e-59         1
## ENSG00000231389    HLA-DPA1         4 3.20211e-53 1.77493e-49         1
## ENSG00000198771       RCSD1         5 8.91780e-43 8.23856e-40         5
## ENSG00000276043       UHRF1         6 1.16106e-45 1.37909e-42         1
## ENSG00000196549         MME         6 4.77191e-49 8.81690e-46         2
##                 wilcox.Top binom.Top    t.p.value wilcox.p.value binom.p.value
##                  <integer> <integer>    <numeric>      <numeric>     <numeric>
## ENSG00000124766          3         1  1.92334e-90    5.72384e-56   1.49550e-41
## ENSG00000128218          4         1 4.45829e-183    1.59904e-79   6.04516e-64
## ENSG00000231389          3         4 5.19868e-214    5.18703e-90   3.20211e-53
## ENSG00000198771          2         1  2.47346e-90    1.43745e-67   8.91780e-43
## ENSG00000276043          6         1 1.76876e-139    6.69601e-56   1.16106e-45
## ENSG00000196549          6         2 8.20175e-139    3.72890e-58   4.77191e-49
```

```r
rm(sce)
```

In addition, multiMarkerStats() will compute a number of new statistics by combining the per-regime statistics. The combined Top value is obtained by simply taking the largest Top value across all tests for a given gene, while the reported p.value is obtained by taking the largest p-value. Ranking on either metric focuses on genes with robust differences that are highly ranked and detected by each of the individual testing regimes. Of course, this might be considered an overly conservative approach in practice, so it is entirely permissible to re-rank the DataFrame according to the Top or p.value for an individual regime (effectively limiting the use of the other regimes’ statistics to diagnostics only).

Write list to file:


```r
tmpFn <- sprintf("%s/%s/Robjects/%s_sce_nz_postDeconv%s%s_%s_clustMarkCombi.Rds",
                 projDir, outDirBit, setName, setSuf, dsiSuf, splSetToGet2)
print(tmpFn)
```

```
## [1] "/ssd/personal/baller01/20200511_FernandesM_ME_crukBiSs2020/AnaWiSce/AnaKmWiC/Robjects/caron_sce_nz_postDeconv_5hCellPerSpl_dsi_PBMMC_ETV6-RUNX1_clustMarkCombi.Rds"
```

```r
saveRDS(combined, file=tmpFn)
rm(combined,tmpFn)
```

## Invalidity of p-values

### From data snooping

All of our DE strategies for detecting marker genes between clusters are statistically flawed to some extent. The DE analysis is performed on the same data used to obtain the clusters, which represents “data dredging” (also known as fishing or data snooping). The hypothesis of interest - are there differences between clusters? - is formulated from the data, so we are more likely to get a positive result when we re-use the data set to test that hypothesis.

The practical effect of data dredging is best illustrated with a simple simulation. We simulate i.i.d. normal values, perform k-means clustering and test for DE between clusters of cells with findMarkers(). The resulting distribution of p-values is heavily skewed towards low values. Thus, we can detect “significant” differences between clusters even in the absence of any real substructure in the data. This effect arises from the fact that clustering, by definition, yields groups of cells that are separated in expression space. Testing for DE genes between clusters will inevitably yield some significant results as that is how the clusters were defined.

Distribution of $p$-values from a DE analysis between two clusters in a simulation with no true subpopulation structure:


```r
library(scran)
set.seed(0)
y <- matrix(rnorm(100000), ncol=200)
clusters <- kmeans(t(y), centers=2)$cluster
out <- findMarkers(y, clusters)
hist(out[[1]]$p.value, col="grey80", xlab="p-value")
```

<img src="clusterMarkerGenes_files/figure-html/unnamed-chunk-53-1.png" width="672" />

For marker gene detection, this effect is largely harmless as the p-values are used only for ranking. However, it becomes an issue when the p-values are used to define “significant differences” between clusters with respect to an error rate threshold. Meaningful interpretation of error rates require consideration of the long-run behavior, i.e., the rate of incorrect rejections if the experiment were repeated many times. The concept of statistical significance for differences between clusters is not applicable if clusters and their interpretations are not stably reproducible across (hypothetical) replicate experiments.

### Nature of replication

The naive application of DE analysis methods will treat counts from the same cluster of cells as replicate observations. This is not the most relevant level of replication when cells are derived from the same biological sample (i.e., cell culture, animal or patient). DE analyses that treat cells as replicates fail to properly model the sample-to-sample variability (Lun and Marioni 2017). The latter is arguably the more important level of replication as different samples will necessarily be generated if the experiment is to be replicated. Indeed, the use of cells as replicates only masks the fact that the sample size is actually one in an experiment involving a single biological sample. This reinforces the inappropriateness of using the marker gene p-values to perform statistical inference.

"We strongly recommend selecting some markers for use in validation studies with an independent replicate population of cells. A typical strategy is to identify a corresponding subset of cells that express the upregulated markers and do not express the downregulated markers. Ideally, a different technique for quantifying expression would also be used during validation, e.g., fluorescent in situ hybridisation or quantitative PCR. This confirms that the subpopulation genuinely exists and is not an artifact of the scRNA-seq protocol or the computational analysis."

See the OSCA chapter on [Marker gene detection](https://osca.bioconductor.org/clustering.html)

**Challenge** Identify markers for a different cluster and try to identify the cell type.

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
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] glue_1.4.2                  Cairo_1.5-12.2             
##  [3] pheatmap_1.0.12             RColorBrewer_1.1-2         
##  [5] dplyr_1.0.5                 scran_1.18.7               
##  [7] scater_1.18.6               SingleCellExperiment_1.12.0
##  [9] SummarizedExperiment_1.20.0 Biobase_2.50.0             
## [11] GenomicRanges_1.42.0        GenomeInfoDb_1.26.7        
## [13] IRanges_2.24.1              S4Vectors_0.28.1           
## [15] BiocGenerics_0.36.1         MatrixGenerics_1.2.1       
## [17] matrixStats_0.58.0          ggplot2_3.3.3              
## [19] knitr_1.32                 
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-7              tools_4.0.3              
##  [3] bslib_0.2.4               utf8_1.2.1               
##  [5] R6_2.5.0                  irlba_2.3.3              
##  [7] vipor_0.4.5               DBI_1.1.1                
##  [9] colorspace_2.0-0          withr_2.4.2              
## [11] tidyselect_1.1.1          gridExtra_2.3            
## [13] compiler_4.0.3            BiocNeighbors_1.8.2      
## [15] DelayedArray_0.16.3       labeling_0.4.2           
## [17] bookdown_0.22             sass_0.3.1               
## [19] scales_1.1.1              stringr_1.4.0            
## [21] digest_0.6.27             rmarkdown_2.7            
## [23] XVector_0.30.0            pkgconfig_2.0.3          
## [25] htmltools_0.5.1.1         sparseMatrixStats_1.2.1  
## [27] limma_3.46.0              highr_0.9                
## [29] rlang_0.4.10              DelayedMatrixStats_1.12.3
## [31] jquerylib_0.1.3           generics_0.1.0           
## [33] farver_2.1.0              jsonlite_1.7.2           
## [35] BiocParallel_1.24.1       RCurl_1.98-1.3           
## [37] magrittr_2.0.1            BiocSingular_1.6.0       
## [39] GenomeInfoDbData_1.2.4    scuttle_1.0.4            
## [41] Matrix_1.3-2              Rcpp_1.0.6               
## [43] ggbeeswarm_0.6.0          munsell_0.5.0            
## [45] fansi_0.4.2               viridis_0.6.0            
## [47] lifecycle_1.0.0           stringi_1.5.3            
## [49] yaml_2.2.1                edgeR_3.32.1             
## [51] zlibbioc_1.36.0           grid_4.0.3               
## [53] dqrng_0.3.0               crayon_1.4.1             
## [55] lattice_0.20-44           cowplot_1.1.1            
## [57] beachmat_2.6.4            locfit_1.5-9.4           
## [59] pillar_1.6.0              igraph_1.2.6             
## [61] codetools_0.2-18          evaluate_0.14            
## [63] vctrs_0.3.7               gtable_0.3.0             
## [65] purrr_0.3.4               assertthat_0.2.1         
## [67] xfun_0.22                 rsvd_1.0.5               
## [69] viridisLite_0.4.0         tibble_3.1.1             
## [71] beeswarm_0.3.1            bluster_1.0.0            
## [73] statmod_1.4.35            ellipsis_0.3.2
```
