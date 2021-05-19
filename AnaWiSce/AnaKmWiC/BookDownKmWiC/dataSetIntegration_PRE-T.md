---
title: "CRUK CI Summer School 2020 - introduction to single-cell RNA-seq analysis"
subtitle: 'Data integration'
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
  dirRel: ".."
  inpDirBit: "AnaWiSce/Ana1"
  outDirBit: "AnaWiSce/Ana1"
  bookType: "mk"
  cacheBool: FALSE  
  setName: "caron"
  splSetToGet: "PRE-T"
  setSuf: "_allCells"
---

# Data integration - PRE-T {#DataIntegrationPretTop}


```r
projDir <- params$projDir
dirRel <- params$dirRel
outDirBit <- params$outDirBit
cacheBool <- params$cacheBool
setName <- params$setName
splSetToGet <- params$splSetToGet
setSuf <- params$setSuf
dsiSuf <- params$dsiSuf # 'dsi' for data set integration

if(params$bookType == "mk"){
	setName <- "caron"
	splSetToGet <- "PRE-T"
	setSuf <- "_allCells"
}
nbPcToComp <- 50
figSize <- 7
```





Source: ['Integrating Datasets'](https://osca.bioconductor.org/integrating-datasets.html) chapter in the OSCA book. Its text is reproduced below with few modifications to adapt it to the data set under scrutiny here.

## Abbreviations
 
* HVG: highly variable genes
* MNN: mutual nearest neighbors
* PBMMC: peripheral blood mononuclear cell
* SCE: SingleCellExperiment

## Motivation

Large single-cell RNA sequencing (scRNA-seq) projects usually need to generate data across multiple batches due to logistical constraints. However, the processing of different batches is often subject to uncontrollable differences, e.g., changes in operator, differences in reagent quality. This results in systematic differences in the observed expression in cells from different batches, which we refer to as “batch effects”. Batch effects are problematic as they can be major drivers of heterogeneity in the data, masking the relevant biological differences and complicating interpretation of the results.

## Load data

Computational correction of these effects is critical for eliminating batch-to-batch variation, allowing data across multiple batches to be combined for common downstream analysis. However, existing methods based on linear models (Ritchie et al. 2015; Leek et al. 2012) assume that the composition of cell populations are either known or the same across batches. To overcome these limitations, bespoke methods have been developed for batch correction of single-cell data (Haghverdi et al. 2018; Butler et al. 2018; Lin et al. 2019) that do not require a priori knowledge about the composition of the population. This allows them to be used in workflows for exploratory analyses of scRNA-seq data where such knowledge is usually unavailable.

## Loading the data

We will load the R file keeping the SCE object with the normalised counts, and subset 1000 cells per sample.


```r
setName <- "caron"
#setSuf <- ""
setSuf <- "_allCells"
tmpFn <- sprintf("%s/%s/Robjects/%s_sce_nz_postDeconv%s.Rds", projDir, outDirBit, setName, setSuf)

print(tmpFn)
```

```
## [1] "/ssd/personal/baller01/20200511_FernandesM_ME_crukBiSs2020/AnaWiSce/AnaKmWiC/Robjects/caron_sce_nz_postDeconv_allCells.Rds"
```

```r
if(!file.exists(tmpFn))
{
	knitr::knit_exit()
}
sce <- readRDS(tmpFn)
sce
```

```
## class: SingleCellExperiment 
## dim: 18431 47830 
## metadata(0):
## assays(2): counts logcounts
## rownames(18431): ENSG00000238009 ENSG00000237491 ... ENSG00000275063
##   ENSG00000271254
## rowData names(11): ensembl_gene_id external_gene_name ... detected
##   gene_sparsity
## colnames: NULL
## colData names(17): Sample Barcode ... cell_sparsity sizeFactor
## reducedDimNames(0):
## altExpNames(0):
```

```r
colnames(rowData(sce))[colnames(rowData(sce)) == "strand"] <- "strandNum"
```

We next subset the data for the PRE-T sample group:


```r
# CaronBourque2020
cb_sampleSheetFn <- file.path(projDir, "Data/CaronBourque2020/SraRunTable.txt")
cb_sampleSheet <- read.table(cb_sampleSheetFn, header=T, sep=",")
splVec <- cb_sampleSheet %>% filter(source_name == splSetToGet) %>%
	pull(Sample.Name) %>% unique

sourceNames <- unique(colData(sce)$source_name)
sceOrig <- sce
sce <- sceOrig[,sce$source_name == splSetToGet ]
nbCells <- 1000
all.sce <- list()
for(spx in splVec)
{
	vec.bc <- colData(sce) %>%
		data.frame() %>%
		filter(Sample.Name == spx) %>%
		sample_n(nbCells) %>%
		pull(Barcode)
	tmpInd <- which(colData(sce)$Barcode %in% vec.bc)
	
	all.sce[[spx]] <- sce[,tmpInd]
}
nbSpl <- length(all.sce)
```

We then apply the standard workflow to each sample separately:
* normalisation,
* variance modelling
* dimensionality reduction
* clustering


```r
#--- normalization ---#
# use logNormCounts()
all.sce <- lapply(all.sce, logNormCounts)

#--- variance-modelling ---#
# model varaince with modelGeneVar()
# find highly variable genes (HVGs) with getTopHVGs()
all.dec <- lapply(all.sce, modelGeneVar)
all.hvgs <- lapply(all.dec, getTopHVGs, prop=0.1)

#--- dimensionality-reduction ---#
# use runPCA()
# then compute embeddings with runTSNE() and runUMAP()
library(BiocSingular)
set.seed(10000)
all.sce <- mapply(FUN=runPCA, x=all.sce, subset_row=all.hvgs, 
    MoreArgs=list(ncomponents=25, BSPARAM=RandomParam()), 
    SIMPLIFY=FALSE)

set.seed(100000)
all.sce <- lapply(all.sce, runTSNE, dimred="PCA")

set.seed(1000000)
all.sce <- lapply(all.sce, runUMAP, dimred="PCA")

#--- clustering ---#
# cluster each sample separately
for (n in names(all.sce)) {
    g <- buildSNNGraph(all.sce[[n]], k=10, use.dimred='PCA')
    clust <- igraph::cluster_walktrap(g)$membership
    #colLabels(all.sce[[n]])  <- factor(clust)
    all.sce[[n]]$label  <- factor(clust)
}
```

To prepare for the batch correction:

* We subset all batches to the common “universe” of features. In this case, it is straightforward as both batches use Ensembl gene annotation.


```r
allNames <- unlist(lapply(all.sce, function(x){rownames(x)}))
allNamesNb <- table(allNames)
universe <- names(allNamesNb)[allNamesNb==nbSpl]
#length(universe)
```

* The size of this common “universe” of features here is the number of features shared by all 2 samples is: 18431.


```r
# Subsetting the SingleCellExperiment object.
uni.sce <- lapply(all.sce, function(x){x[universe,]})
# Also subsetting the variance modelling results, for convenience.
uni.dec <- lapply(all.dec, function(x){x[universe,]})
```

* We rescale each batch to adjust for differences in sequencing depth between batches. The multiBatchNorm() function recomputes log-normalized expression values after adjusting the size factors for systematic differences in coverage between SingleCellExperiment (SCE) objects. (Size factors only remove biases between cells within a single batch.) This improves the quality of the correction by removing one aspect of the technical differences between batches.


```r
# rescale each batch to adjust for differences in sequencing depth between batches
rescaled <- multiBatchNorm(uni.sce, batch = "Sample.Name")
```

* We perform feature selection by averaging the variance components across all batches with the combineVar() function. We compute the average as it is responsive to batch-specific HVGs while still preserving the within-batch ranking of genes.


```r
# compute average variance components across samples
#combined.dec <- combineVar(uni.dec[[1]], uni.dec[[2]], uni.dec[[3]], uni.dec[[4]])
combined.dec <- combineVar(uni.dec)
# identify highly variables genes
# here as those with a positive biological component
chosen.hvgs <- combined.dec$bio > 0
#sum(chosen.hvgs)
```

Number of HVGs: 10485.

When integrating datasets of variable composition, it is generally safer to err on the side of including more genes than are used in a single dataset analysis, to ensure that markers are retained for any dataset-specific subpopulations that might be present. For a top X selection, this means using a larger X (say, ~5000), or in this case, we simply take all genes above the trend.

Alternatively, a more forceful approach to feature selection can be used based on marker genes from within-batch comparisons.

## Diagnosing batch effects

Before we actually perform any correction, it is worth examining whether there is any batch effect in this dataset. We combine the SCE objects and perform a PCA on the log-expression values for all genes with positive (average) biological components.


```r
# Synchronizing the metadata for cbind()ing.
for (i in 2:nbSpl)
{
  identical(rowData(rescaled[[1]]), rowData(rescaled[[i]]))
}

rescaled[[1]]$batch <- rescaled[[1]]$Sample.Name
rescaled2 <- lapply(rescaled, function(x){x$batch <- x$Sample.Name; x})
rescaled <- rescaled2

# concat matrices:
uncorrected <- do.call(cbind, rescaled)

# Perform PCA
# Using RandomParam() as it is more efficient for file-backed matrices.
set.seed(0010101010)
uncorrected <- runPCA(uncorrected,
                      subset_row=chosen.hvgs,
                      BSPARAM=BiocSingular::RandomParam())
```

We use graph-based clustering on the components to obtain a summary of the population structure.

As the samples should be replicates, each cluster should ideally consist of cells from each batch. However, we instead see clusters that are comprised of cells from a single batch. This indicates that cells of the same type are artificially separated due to technical differences between batches.


```r
# build shared nearest-neighbour graph
snn.gr <- buildSNNGraph(uncorrected, use.dimred="PCA")
# identify cluster with the walk trap method
clusters <- igraph::cluster_walktrap(snn.gr)$membership
# get number of cells for each {cluster, batch} pair
tab <- table(Cluster=clusters, Batch=uncorrected$batch)
#tab
```

```r
tmpMat <- data.frame("clusters"=clusters, "batch"=uncorrected$batch)
```



Cluster size and cell contribution by sample:


```r
tmpMatTab <- table(tmpMat)
sortVecNames <- tmpMatTab %>% rowSums %>% sort(decreasing=TRUE) %>% names
tmpMat$clusters <- factor(tmpMat$clusters, levels=sortVecNames)
tmpMatDf <- tmpMatTab[sortVecNames,] %>% data.frame()
p1 <- ggplot(data=tmpMatDf, aes(x=clusters,y=Freq, fill=batch)) +
	geom_col() +
  ggtitle("uncorrected, cell numbers")
p2 <- ggplot(data=tmpMat, aes(x=clusters, fill=batch)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  ggtitle("uncorrected, proportions")
```

```r
gridExtra::grid.arrange(p1, p2)
```

<img src="dataSetIntegration_PRE-T_files/figure-html/unnamed-chunk-15-1.png" width="672" />

We can also visualize the uncorrected coordinates using a t-SNE plot. The strong separation between cells from different batches is consistent with the clustering results.


```r
set.seed(1111001)
uncorrected <- runTSNE(uncorrected, dimred="PCA")
plotTSNE(uncorrected, colour_by="batch")
```

<img src="dataSetIntegration_PRE-T_files/figure-html/unnamed-chunk-16-1.png" width="672" />

Of course, the other explanation for batch-specific clusters is that there are cell types that are unique to each batch. The degree of intermingling of cells from different batches is not an effective diagnostic when the batches involved might actually contain unique cell subpopulations. If a cluster only contains cells from a single batch, one can always debate whether that is caused by a failure of the correction method or if there is truly a batch-specific subpopulation. For example, do batch-specific metabolic or differentiation states represent distinct subpopulations? Or should they be merged together? We will not attempt to answer this here, only noting that each batch correction algorithm will make different (and possibly inappropriate) decisions on what constitutes “shared” and “unique” populations.

## Linear regression

Batch effects in bulk RNA sequencing studies are commonly removed with linear regression. This involves fitting a linear model to each gene’s expression profile, setting the undesirable batch term to zero and recomputing the observations sans the batch effect, yielding a set of corrected expression values for downstream analyses. Linear modelling is the basis of the `removeBatchEffect()` function from the limma package (Ritchie et al. 2015) as well the `comBat()` function from the sva package (Leek et al. 2012).

To use this approach in a scRNA-seq context, we assume that the composition of cell subpopulations is the same across batches. We also assume that the batch effect is additive, i.e., any batch-induced fold-change in expression is the same across different cell subpopulations for any given gene. These are strong assumptions as batches derived from different individuals will naturally exhibit variation in cell type abundances and expression. Nonetheless, they may be acceptable when dealing with batches that are technical replicates generated from the same population of cells. (In fact, when its assumptions hold, linear regression is the most statistically efficient as it uses information from all cells to compute the common batch vector.) Linear modelling can also accommodate situations where the composition is known a priori by including the cell type as a factor in the linear model, but this situation is even less common.

We use the `rescaleBatches()` function from the `batchelor` package to remove the batch effect. This is roughly equivalent to applying a linear regression to the log-expression values per gene, with some adjustments to improve performance and efficiency. For each gene, the mean expression in each batch is scaled down until it is equal to the lowest mean across all batches. We deliberately choose to scale all expression values down as this mitigates differences in variance when batches lie at different positions on the mean-variance trend. (Specifically, the shrinkage effect of the pseudo-count is greater for smaller counts, suppressing any differences in variance across batches.) An additional feature of `rescaleBatches()` is that it will preserve sparsity in the input matrix for greater efficiency, whereas other methods like `removeBatchEffect()` will always return a dense matrix.


```r
rescaled2 <- rescaleBatches(rescaled)
rescaled2
```

```
## class: SingleCellExperiment 
## dim: 18431 2000 
## metadata(0):
## assays(1): corrected
## rownames(18431): ENSG00000000003 ENSG00000000419 ... ENSG00000285486
##   ENSG00000285492
## rowData names(0):
## colnames: NULL
## colData names(1): batch
## reducedDimNames(0):
## altExpNames(0):
```

After clustering, we observe that most clusters consist of mixtures of cells from the two replicate batches, consistent with the removal of the batch effect. This conclusion is supported by the apparent mixing of cells from different batches in Figure 13.2. However, at least one batch-specific cluster is still present, indicating that the correction is not entirely complete. This is attributable to violation of one of the aforementioned assumptions, even in this simple case involving replicated batches.


```r
set.seed(1010101010) # To ensure reproducibility of IRLBA.
rescaled2 <- runPCA(rescaled2, subset_row=chosen.hvgs, exprs_values="corrected")

snn.gr <- buildSNNGraph(rescaled2, use.dimred="PCA")
clusters.resc <- igraph::cluster_walktrap(snn.gr)$membership
tab.resc <- table(Cluster=clusters.resc, Batch=rescaled2$batch)
#tab.resc
```


```r
tmpMat <- data.frame("clusters"=clusters.resc, "batch"=rescaled2$batch)
```

Cluster size and cell contribution by sample, with clusters sorted by size:


```r
tmpMatTab <- table(tmpMat)
sortVecNames <- tmpMatTab %>% rowSums %>% sort(decreasing=TRUE) %>% names
tmpMat$clusters <- factor(tmpMat$clusters, levels=sortVecNames)
tmpMatDf <- tmpMatTab[sortVecNames,] %>% data.frame()
p1 <- ggplot(data=tmpMatDf, aes(x=clusters,y=Freq, fill=batch)) +
	geom_col()
p2 <- ggplot(data=tmpMat, aes(x=clusters, fill=batch)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent)
```

```r
gridExtra::grid.arrange(p1, p2)
```

<img src="dataSetIntegration_PRE-T_files/figure-html/unnamed-chunk-21-1.png" width="672" />

Compute and plot t-SNE:


```r
rescaled2 <- runTSNE(rescaled2, dimred="PCA")
rescaled2$batch <- factor(rescaled2$batch)
plotTSNE(rescaled2, colour_by="batch")
```

<img src="dataSetIntegration_PRE-T_files/figure-html/unnamed-chunk-22-1.png" width="672" />



## Performing MNN correction

### Algorithm overview

Consider a cell a in batch A, and identify the cells in batch B that are nearest neighbors to a in the expression space defined by the selected features. Repeat this for a cell b in batch B, identifying its nearest neighbors in A. Mutual nearest neighbors are pairs of cells from different batches that belong in each other’s set of nearest neighbors. The reasoning is that MNN pairs represent cells from the same biological state prior to the application of a batch effect - see Haghverdi et al. (2018) for full theoretical details. Thus, the difference between cells in MNN pairs can be used as an estimate of the batch effect, the subtraction of which yields batch-corrected values.

Compared to linear regression, MNN correction does not assume that the population composition is the same or known beforehand. This is because it learns the shared population structure via identification of MNN pairs and uses this information to obtain an appropriate estimate of the batch effect. Instead, the key assumption of MNN-based approaches is that the batch effect is orthogonal to the biology in high-dimensional expression space. Violations reduce the effectiveness and accuracy of the correction, with the most common case arising from variations in the direction of the batch effect between clusters. Nonetheless, the assumption is usually reasonable as a random vector is very likely to be orthogonal in high-dimensional space.


### Application to the data

The `batchelor` package provides an implementation of the MNN approach via the `fastMNN()` function. (Unlike the MNN method originally described by Haghverdi et al. (2018), the `fastMNN()` function performs PCA to reduce the dimensions beforehand and speed up the downstream neighbor detection steps.) We apply it to our two PBMC batches to remove the batch effect across the highly variable genes in `chosen.hvgs`. To reduce computational work and technical noise, all cells in all batches are projected into the low-dimensional space defined by the top d principal components. Identification of MNNs and calculation of correction vectors are then performed in this low-dimensional space.


```r
# Using randomized SVD here, as this is faster than 
# irlba for file-backed matrices.
set.seed(1000101001)
mnn.out <- fastMNN(rescaled,
                   auto.merge=TRUE,
                   d=50,
                   k=20,
                   subset.row=chosen.hvgs,
                   BSPARAM=BiocSingular::RandomParam(deferred=TRUE))
mnn.out
```

```
## class: SingleCellExperiment 
## dim: 10485 2000 
## metadata(2): merge.info pca.info
## assays(1): reconstructed
## rownames(10485): ENSG00000000003 ENSG00000000457 ... ENSG00000285444
##   ENSG00000285476
## rowData names(1): rotation
## colnames: NULL
## colData names(1): batch
## reducedDimNames(1): corrected
## altExpNames(0):
```

```r
mnn.out.corre.dim <- dim(reducedDim(mnn.out, "corrected"))
mnn.out.corre.dim
```

```
## [1] 2000   50
```

```r
mnn.out.recon.dim <- dim(assay(mnn.out, "reconstructed"))
mnn.out.recon.dim
```

```
## [1] 10485  2000
```

The function returns a SCE object containing corrected values for downstream analyses like clustering or visualization. Each column of `mnn.out` corresponds to a cell in one of the batches, while each row corresponds to an input gene in `chosen.hvgs`. The batch field in the column metadata contains a vector specifying the batch of origin of each cell.

The `corrected` matrix in the `reducedDims()` contains the low-dimensional corrected coordinates for all cells, which we will use in place of the PCs in our downstream analyses (2000 cells and 50 PCs).

A `reconstructed` matrix in the `assays()` contains the corrected expression values for each gene in each cell, obtained by projecting the low-dimensional coordinates in corrected back into gene expression space (10485 genes and 2000 cells). We do not recommend using this for anything other than visualization.


```r
print(assay(mnn.out, "reconstructed")[1:5,1:3])
```

```
## <5 x 3> matrix of class LowRankMatrix and type "double":
##                          [,1]          [,2]          [,3]
## ENSG00000000003  5.679825e-04  5.770801e-04 -4.069776e-04
## ENSG00000000457  8.492035e-04 -2.059444e-04  1.583788e-04
## ENSG00000000938 -1.399295e-03 -1.799694e-03 -2.562064e-03
## ENSG00000000971 -3.513720e-05 -4.466999e-05 -9.200578e-05
## ENSG00000001036 -6.472921e-05 -1.734277e-04 -6.579931e-04
```

The most relevant parameter for tuning `fastMNN()` is `k`, which specifies the number of nearest neighbors to consider when defining MNN pairs. This can be interpreted as the minimum anticipated frequency of any shared cell type or state in each batch. Increasing `k` will generally result in more aggressive merging as the algorithm is more generous in matching subpopulations across batches. It can occasionally be desirable to increase `k` if one clearly sees that the same cell types are not being adequately merged across batches.

<!--
See Chapter 32 for an example of a more complex fastMNN() merge involving several human pancreas datasets generated by different authors on different patients with different technologies.
-->

## Correction diagnostics

### Mixing between batches

We cluster on the low-dimensional corrected coordinates to obtain a partitioning of the cells that serves as a proxy for the population structure. If the batch effect is successfully corrected, clusters corresponding to shared cell types or states should contain cells from multiple batches. We see that all clusters contain contributions from each batch after correction, consistent with our expectation that the batches are replicates of each other.


```r
snn.gr <- buildSNNGraph(mnn.out, use.dimred="corrected", k=20)
clusters.mnn <- igraph::cluster_walktrap(snn.gr)$membership
tab.mnn <- table(Cluster=clusters.mnn, Batch=mnn.out$batch)
tab.mnn
```

```
##        Batch
## Cluster GSM3872440 GSM3872441
##      1         183         37
##      2           2         74
##      3           1         30
##      4         438        241
##      5         116         31
##      6          15        322
##      7           7        107
##      8           8         43
##      9           2         58
##      10        228         57
```

Cluster size and cell contribution by sample, with clusters sorted by size:


```r
tmpMat <- data.frame("clusters"=clusters.mnn, "batch"=mnn.out$batch)
tmpMatTab <- table(tmpMat)
sortVecNames <- tmpMatTab %>% rowSums %>% sort(decreasing=TRUE) %>% names
tmpMat$clusters <- factor(tmpMat$clusters, levels=sortVecNames)
tmpMatTab <- table(tmpMat)
tmpMatDf <- tmpMatTab[sortVecNames,] %>% data.frame()
p1 <- ggplot(data=tmpMatDf, aes(x=clusters,y=Freq, fill=batch)) +
	geom_col()
p2 <- ggplot(data=tmpMat, aes(x=clusters, fill=batch)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent)
```

```r
gridExtra::grid.arrange(p1, p2)
```

<img src="dataSetIntegration_PRE-T_files/figure-html/unnamed-chunk-27-1.png" width="672" />

We can also compute the variation in the log-abundances to rank the clusters with the greatest variability in their proportional abundances across batches. We can then focus on batch-specific clusters that may be indicative of incomplete batch correction. Obviously, though, this diagnostic is subject to interpretation as the same outcome can be caused by batch-specific populations; some prior knowledge about the biological context is necessary to distinguish between these two possibilities. The table below shows the number of cells for each cluster (row) and sample (column) together with the variance in cell number across these samples ('var' column). 

Also bear in mind that the variance is computed across 2 samples here and only serves to sort clusters.


```r
# Avoid minor difficulties with the 'table' class.
tab.mnn <- unclass(tab.mnn)

# Using a large pseudo.count to avoid unnecessarily
# large variances when the counts are low.
norm <- normalizeCounts(tab.mnn, pseudo_count=10)

# Ranking clusters by the largest variances.
rv <- rowVars(norm) %>% round(2)

# show
#DataFrame(Batch=tab.mnn, var=rv)[order(rv, decreasing=TRUE),]
DataFrame(tab.mnn, var=rv)[order(rv, decreasing=TRUE),]
```

```
## DataFrame with 10 rows and 3 columns
##    GSM3872440 GSM3872441       var
##     <integer>  <integer> <numeric>
## 6          15        322      6.96
## 2           2         74      3.94
## 7           7        107      3.87
## 9           2         58      3.13
## 1         183         37      2.08
## 3           1         30      1.73
## 10        228         57      1.67
## 5         116         31      1.31
## 8           8         43      1.21
## 4         438        241      0.35
```

We can also visualize the corrected coordinates using a t-SNE plot. The presence of visual clusters containing cells from both batches provides a comforting illusion that the correction was successful.


```r
set.seed(0010101010)
mnn.out <- runTSNE(mnn.out, dimred="corrected")

mnn.out$batch <- factor(mnn.out$batch)
plotTSNE(mnn.out, colour_by="batch")
```

<img src="dataSetIntegration_PRE-T_files/figure-html/unnamed-chunk-29-1.png" width="672" />

```r
#mnn.out$type <- gsub("_[1-4]","",mnn.out$batch)
#p <- plotTSNE(mnn.out, colour_by="batch", shape_by="type")
#p + facet_wrap(. ~ mnn.out$type)
```

For `fastMNN()`, one useful diagnostic is the proportion of variance within each batch that is lost during MNN correction. Specifically, this refers to the within-batch variance that is removed during orthogonalization with respect to the average correction vector at each merge step. This is returned via the `lost.var` field in the metadata of `mnn.out`, which contains a matrix of the variance lost in each batch (column) at each merge step (row).


```r
round(metadata(mnn.out)$merge.info$lost.var,2)
```

```
##      GSM3872440 GSM3872441
## [1,]       0.01       0.01
```
Large proportions of lost variance (>10%) suggest that correction is removing genuine biological heterogeneity. This would occur due to violations of the assumption of orthogonality between the batch effect and the biological subspace (Haghverdi et al. 2018). In this case, the proportion of lost variance is small, indicating that non-orthogonality is not a major concern.

The following t-SNE shows the clusters identified:


```r
mnn.out$cluster <- paste0("c", clusters.mnn)
p <- plotTSNE(mnn.out, colour_by="cluster")
p
```

<img src="dataSetIntegration_PRE-T_files/figure-html/unnamed-chunk-31-1.png" width="672" />


```r
p + facet_wrap(~colData(mnn.out)$batch)
```

<img src="dataSetIntegration_PRE-T_files/figure-html/unnamed-chunk-32-1.png" width="672" />

The following t-SNE plots show expression levels of known cell type marker genes.


```r
genesToShow <- c(
		 "CD79A", # CD79A 	B ***
		 "CST3", # CST3 	monocytes ***
		 "CD3D", # CD3D 	 T cells ***
		 "HBA1" # HBA1 	 erythrocytes ***
	 	)

tmpInd <- which(rowData(uncorrected)$Symbol %in% genesToShow)
ensToShow <- rowData(uncorrected)$ensembl_gene_id[tmpInd]
```

B cells:


```r
genex <- ensToShow[1]
	p <- plotTSNE(mnn.out, colour_by = genex, by_exprs_values="reconstructed")
	p <- p + ggtitle(
			paste("B cells", genex,
			rowData(uncorrected)[genex,"Symbol"])
		)
	print(p)
```

<img src="dataSetIntegration_PRE-T_files/figure-html/unnamed-chunk-34-1.png" width="672" />

T cells:


```r
genex <- ensToShow[3]
	p <- plotTSNE(mnn.out, colour_by = genex, by_exprs_values="reconstructed")
	p <- p + ggtitle(
			paste("T cells", genex,
			rowData(uncorrected)[genex,"Symbol"])
		)
	print(p)
```

<img src="dataSetIntegration_PRE-T_files/figure-html/unnamed-chunk-35-1.png" width="672" />

monocytes:


```r
genex <- ensToShow[2]
	p <- plotTSNE(mnn.out, colour_by = genex, by_exprs_values="reconstructed")
	p <- p + ggtitle(
			paste("monocytes", genex,
			rowData(uncorrected)[genex,"Symbol"])
		)
	print(p)
```

<img src="dataSetIntegration_PRE-T_files/figure-html/unnamed-chunk-36-1.png" width="672" />

erythrocytes:


```r
genex <- ensToShow[4]
	p <- plotTSNE(mnn.out, colour_by = genex, by_exprs_values="reconstructed")
	p <- p + ggtitle(
			paste("erythrocytes", genex,
			rowData(uncorrected)[genex,"Symbol"])
		)
	print(p)
```

<img src="dataSetIntegration_PRE-T_files/figure-html/unnamed-chunk-37-1.png" width="672" />

Compare to the uncorrected values, T cells:


```r
genex <- ensToShow[3]
	p <- plotTSNE(uncorrected, colour_by = genex)
	p <- p + ggtitle(
			paste("T cells", genex,
			rowData(uncorrected)[genex,"Symbol"])
		)
	print(p)
```

<img src="dataSetIntegration_PRE-T_files/figure-html/unnamed-chunk-38-1.png" width="672" />
Compare to the uncorrected values, erythrocytes:


```r
genex <- ensToShow[4]
	p <- plotTSNE(uncorrected, colour_by = genex)
	p <- p + ggtitle(
			paste("erythrocytes", genex,
			rowData(uncorrected)[genex,"Symbol"])
		)
	print(p)
```

<img src="dataSetIntegration_PRE-T_files/figure-html/unnamed-chunk-39-1.png" width="672" />

Other genes (exercise)


```r
genesToShow2 <- c(
		 "IL7R", # IL7R, CCR7 	Naive CD4+ T
		 "CCR7", # IL7R, CCR7 	Naive CD4+ T
		 "S100A4", # IL7R, S100A4 	Memory CD4+
		 "CD14", # CD14, LYZ 	CD14+ Mono
		 "LYZ", # CD14, LYZ 	CD14+ Mono
		 "MS4A1", # MS4A1 	B
		 "CD8A", # CD8A 	CD8+ T
		 "FCGR3A", # FCGR3A, MS4A7 	FCGR3A+ Mono
		 "MS4A7", # FCGR3A, MS4A7 	FCGR3A+ Mono
		 "GNLY", # GNLY, NKG7 	NK
		 "NKG7", # GNLY, NKG7 	NK
		 "FCER1A", # DC
		 "CST3", # DC
		 "PPBP" # Platelet
		)
```


```r
tmpInd <- which(rowData(uncorrected)$Symbol %in% genesToShow2)
ensToShow <- rowData(uncorrected)$ensembl_gene_id[tmpInd]
table(ensToShow %in% rownames(rowData(mnn.out)))
ensToShow <- ensToShow[ensToShow %in% rownames(rowData(mnn.out))]
```



```r
for (genex in ensToShow)
{
	p <- plotTSNE(mnn.out, colour_by = genex, by_exprs_values="reconstructed")
	p <- p + ggtitle(paste(genex, rowData(uncorrected)[genex,"Symbol"]))
	print(p)
}
```

### Preserving biological heterogeneity

#### Comparison between within-batch clusters and across-batch clusters obtained after MNN correction

Another useful diagnostic check is to compare the clustering within each batch to the clustering of the merged data. Accurate data integration should preserve variance within each batch as there should be nothing to remove between cells in the same batch. This check complements the previously mentioned diagnostics that only focus on the removal of differences between batches. Specifically, it protects us against cases where the correction method simply aggregates all cells together, which would achieve perfect mixing but also discard the biological heterogeneity of interest.

Ideally, we should see a many-to-1 mapping where the across-batch clustering is nested inside the within-batch clusterings. This indicates that any within-batch structure was preserved after correction while acknowledging that greater resolution is possible with more cells. In practice, more discrepancies can be expected even when the correction is perfect, due to the existence of closely related clusters that were arbitrarily separated in the within-batch clustering. As a general rule, we can be satisfied with the correction if the vast majority of entries are zero, though this may depend on whether specific clusters of interest are gained or lost.

One heatmap is generated for each dataset, where each entry is colored according to the number of cells with each pair of labels (before and after correction), on the log10 scale with pseudocounts (+10) for a smoother color transition (so a minimum value of log10(0+10) == 1). 










```r
plotList <- vector(mode = "list", length = length(splVec))
treeList <- vector(mode = "list", length = length(splVec))
for (splIdx in 1:length(splVec)) {
  # heatmap
  tab <- table(
    paste("before", colLabels(rescaled[[splIdx]]), sep="_"),
    paste("after", clusters.mnn[rescaled2$batch==splVec[splIdx]], sep="_")
    )
  plotList[[splIdx]] <- pheatmap(log10(tab+10),
                                 cluster_row=FALSE,
                                 cluster_col=FALSE,
                                 col=rev(viridis::magma(100)),
                                 main=sprintf("%s",
                                              splVec[splIdx]),
                                 silent=TRUE,
                                 fontsize=7)
  # cluster tree:
  combined <- cbind(
    cl.1=colLabels(rescaled[[splIdx]]),
    cl.2=clusters.mnn[rescaled2$batch==splVec[splIdx]])
  treeList[[splIdx]]  <- clustree(combined, prefix="cl.", edge_arrow=FALSE) +
    ggtitle(splVec[splIdx]) +
    #theme(legend.background = element_rect(color = "yellow")) +
    #theme(legend.position='bottom') +
    #theme(legend.box="vertical") +
    #theme(legend.box="horizontal") +
    theme(legend.margin=margin()) #+
    #guides(fill=guide_legend(nrow=2, byrow=FALSE))
    #theme(legend.position = "none")
}
```


```r
g_legend<-function(a.gplot){
   tmp <- ggplot_gtable(ggplot_build(a.gplot))
   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
   legend <- tmp$grobs[[leg]]
   return(legend)
}

redrawClutree <- function(p){
#p <- treeList[[1]] + theme(legend.position='bottom')
#p <- p + theme(legend.background = element_rect(color = "yellow"))
p <- p + theme(legend.justification = "left")
#p <- p + theme(legend.justification = c(0,1))
#lemon::gtable_show_names(p)
pNoLeg <- p + theme(legend.position = "none")
# edge colour:
pEdgeCol <- p +
  #guides(edge_colour = FALSE) +
  guides(edge_alpha = FALSE) +
  guides(size = FALSE) +
  guides(colour = FALSE) 
pEdgeCol.leg <- g_legend(pEdgeCol)
# edge alpha:
pEdgeAlpha <- p +
  guides(edge_colour = FALSE) +
  #guides(edge_alpha = FALSE) +
  guides(size = FALSE) +
  guides(colour = FALSE) 
pEdgeAlpha.leg <- g_legend(pEdgeAlpha)
# size
pSize <- p +
  guides(edge_colour = FALSE) +
  guides(edge_alpha = FALSE) +
  #guides(size = FALSE) +
  guides(colour = FALSE) 
pSize.leg <- g_legend(pSize)
# colour
pColour <- p +
  guides(edge_colour = FALSE) +
  guides(edge_alpha = FALSE) +
  guides(size = FALSE) #+
  #guides(colour = FALSE) 
pColour.leg <- g_legend(pColour)

#gridExtra::grid.arrange(pNoLeg, pEdgeCol.leg, nrow=2, ncol=1, heights=c(unit(.8, "npc"), unit(.2, "npc")))
if(FALSE)
{
grobx <- gridExtra::grid.arrange(pNoLeg,
                        pEdgeCol.leg,
                        pEdgeAlpha.leg,
                        pColour.leg,
                        pSize.leg,
                        nrow=3, ncol=2,
                        heights=c(unit(.8, "npc"),
                                  unit(.1, "npc"),
                                  unit(.1, "npc")),
                        widths=c(unit(.3, "npc"), unit(.7, "npc")),
                        layout_matrix=matrix(c(1,1,2,5,4,3), ncol=2, byrow=TRUE)
                        )
}
if(FALSE)
{
grobx <- gridExtra::arrangeGrob(pNoLeg,
                        pEdgeCol.leg,
                        pEdgeAlpha.leg,
                        pColour.leg,
                        pSize.leg,
                        #nrow=3, ncol=2,
                        #layout_matrix=matrix(c(1,1,2,5,4,3), ncol=2, byrow=TRUE),
                        nrow=2, ncol=3,
                        layout_matrix=matrix(c(1,1,2,5,4,3), ncol=3, byrow=FALSE),
                        widths=c(unit(.70, "npc"),
                                  unit(.15, "npc"),
                                  unit(.15, "npc")),
                        heights=c(unit(.7, "npc"),
                                 unit(.3, "npc"))
                        )
}

grobx <- gridExtra::arrangeGrob(pNoLeg,
                        pEdgeCol.leg,
                        pEdgeAlpha.leg,
                        #pColour.leg,
                        pSize.leg,
                        nrow=1, ncol=4,
                        layout_matrix=matrix(c(1,2,3,4), ncol=4, byrow=TRUE),
                        widths=c(unit(.64, "npc"),
                                  unit(.12, "npc"),
                                  unit(.12, "npc"),
                                  unit(.12, "npc"))
                        )
}
##gx <- redrawClutree(treeList[[1]] + theme(legend.position='bottom'))
##grid::grid.draw(gx)
## fine # gxList <- lapply(treeList, function(x){redrawClutree(x+theme(legend.position='bottom'))})
gxList <- lapply(treeList, function(x){redrawClutree(x)})
##gridExtra::marrangeGrob(gxList, nrow=2, ncol=2)
```






```r
grobList <- lapply(plotList, function(x){x[[4]]})
gridExtra::grid.arrange(grobs = grobList,
      ncol=2,
      top = grid::textGrob("clusterings concordance (number of cells, log10 scale)",
                     gp=grid::gpar(fontsize=12,font=3))
)
```

<img src="dataSetIntegration_PRE-T_files/figure-html/unnamed-chunk-50-1.png" width="672" />

The redistribution of cells from one set of clusters to another, here 'within-batch before' and 'across-batch after' correction, may also be visualized with a clustering tree [clustree](https://cran.r-project.org/package=clustree). Clusters are represented as filled circles colored by cluster set ('before' in pink, 'after' in blue) and sized by cell number. A pair of clusters from two sets are linked according to the number of cells they share with a link that informs on the number of cells shared (color) and the 'incoming node' proportion for the node it points to (transparency). Although these plots convey more information than heatmaps below, they may not be as easy to read.




```r
figSize <- 7
```


```r
#```{r, fig.height=figSize*length(treeList)/2, fig.width=figSize}
#gridExtra::grid.arrange(grobs = treeList,
gridExtra::grid.arrange(grobs = gxList,
      ncol=1
)
```

<img src="dataSetIntegration_PRE-T_files/figure-html/unnamed-chunk-53-1.png" width="672" />

The same plots in more compact form with no legend:

<!-- remove legend and have plots in two columns -->


```r
treeList <- lapply(treeList, function(p){
  p +
    guides(edge_colour = FALSE) +
    guides(edge_alpha = FALSE) +
    guides(size = FALSE) +
    guides(colour = FALSE) 
})
```






```r
gridExtra::grid.arrange(grobs = treeList,
      ncol=2
)
```

<img src="dataSetIntegration_PRE-T_files/figure-html/unnamed-chunk-56-1.png" width="672" />



```r
knitr::knit_exit()
```

#### Coassignment probabilities

Another evaluation approach is to compute the **coassignment probabilities**, i.e. the probability that cells from two within-batch clusters are clustered together in the across-batch clustering. High probabilities off the diagonal indicate that within-batch clusters are merged in the across-batch analysis. We would generally expect low off-diagonal probabilities for most pairs of clusters, though this may not be reasonably possible if the within-batch clusters were poorly separated in the first place.

The plots below display the coassignment probabilities for the within-batch clusters, based on coassignment of cells in the across-batch clusters obtained after MNN correction. One heatmap is generated for each sample, where each entry is colored according to the coassignment probability between each pair of within-batch clusters:




```r
# coassignProb manual: now deprecated for pairwiseRand. 
# Note that the coassignment probability is closely related to the Rand index-based ratios broken down by cluster pair in pairwiseRand with mode="ratio" and adjusted=FALSE. The off-diagonal coassignment probabilities are simply 1 minus the off-diagonal ratio, while the on-diagonal values differ only by the lack of consideration of pairs of the same cell in pairwiseRand. 

plotList <- vector(mode = "list", length = length(splVec))
for (splIdx in 1:length(splVec)) {
  tab <- coassignProb(colLabels(rescaled[[splIdx]]),
                      clusters.mnn[rescaled2$batch==splVec[splIdx]])
  plotList[[splIdx]] <- pheatmap(tab,
                                 cluster_row=FALSE,
                                 cluster_col=FALSE,
                                 col=rev(viridis::magma(100)),
                                 main=sprintf("%s probabilities", splVec[splIdx]),
                                 silent=TRUE)
}
grobList <- lapply(plotList, function(x){x[[4]]})
gridExtra::grid.arrange(grobs = grobList,
      ncol=2
)
```

<img src="dataSetIntegration_PRE-T_files/figure-html/unnamed-chunk-59-1.png" width="672" />

Note that the coassignment probability is closely related to the Rand index-based ratios broken down by cluster pair (in `pairwiseRand()` with mode="ratio" and adjusted=FALSE). The Rand index is introduced below.
 



#### Rand index

Finally, we can summarize the agreement between clusterings by computing the **Rand index**. This provides a simple metric that we can use to assess the preservation of variation by different correction methods. Larger rand indices (i.e., closer to 1) are more desirable, though this must be balanced against the ability of each method to actually remove the batch effect.







```r
# pairwiseRand(), index, adjusted
ariVec <- vector(mode = "numeric", length = length(splVec))
names(ariVec) <- splVec
for (splIdx in 1:length(splVec)) {
  ariVec[splIdx] <- pairwiseRand(
    ref=as.integer(clusters.mnn[rescaled2$batch==splVec[splIdx]]),
    alt=as.integer(colLabels(rescaled[[splIdx]])),
    mode="index")
}
ariVec <- round(ariVec,2)
ariVec
```

```
## GSM3872440 GSM3872441 
##       0.64       0.57
```

A sample may show a low Rand index value if cells grouped together in a small cluster before correction are split into distinct clusters after correction because the latter comprise cell populations not observed in that sample but present in other samples.

This would be the case of GSM3872440 with far fewer erythrocytes (group in a single cluster) than GSM3872441, in which subtypes can be distinguished.
<!--
-->



We can also break down the **adjusted Rand index (ARI)** into per-cluster ratios for more detailed diagnostics. For example, we could see low ratios off the diagonal if distinct clusters in the within-batch clustering were incorrectly aggregated in the merged clustering. Conversely, we might see low ratios on the diagonal if the correction inflated or introduced spurious heterogeneity inside a within-batch cluster.




```r
# pairwiseRand(), ratio, adjusted
# square numeric matrix is returned with number of rows equal to the number of unique levels in ref.

tabList <- vector(mode = "list", length = length(splVec))
for (splIdx in 1:length(splVec)) {
  tabList[[splIdx]] <- pairwiseRand(
	  ref=as.integer(colLabels(rescaled[[splIdx]])),
    alt=as.integer(clusters.mnn[rescaled2$batch==splVec[splIdx]])
	)
}
randVal <- unlist(tabList) 

## make breaks from combined range
limits <- c(
  min(randVal, na.rm = TRUE),
  max(randVal, na.rm = TRUE))
limits <- quantile(randVal, probs=c(0.05, 0.95), na.rm = TRUE)

Breaks <- seq(limits[1], limits[2],
              length = 100)

plotList <- vector(mode = "list", length = length(splVec))
for (splIdx in 1:length(splVec)) {
  plotList[[splIdx]] <- pheatmap(tabList[[splIdx]],
                                 cluster_row=FALSE,
                                 cluster_col=FALSE,
                                 col=rev(viridis::magma(100)),
                                 breaks=Breaks,
                                 main=sprintf("%s ratio", splVec[splIdx]),
                                 silent=TRUE)
}
grobList <- lapply(plotList, function(x){x[[4]]})
gridExtra::grid.arrange(grobs = grobList,
      ncol=2
)
```

<img src="dataSetIntegration_PRE-T_files/figure-html/unnamed-chunk-66-1.png" width="672" />

## Encouraging consistency with marker genes

In some situations, we will already have performed within-batch analyses to characterize salient aspects of population heterogeneity. This is not uncommon when merging datasets from different sources where each dataset has already been analyzed, annotated and interpreted separately. It is subsequently desirable for the integration procedure to retain these “known interesting” aspects of each dataset in the merged dataset. We can encourage this outcome by using the marker genes within each dataset as our selected feature set for `fastMNN()` and related methods. This focuses on the relevant heterogeneity and represents a semi-supervised approach that is a natural extension of the strategy described in Section 8.4.

We identify the top marker genes from pairwise Wilcoxon ranked sum tests between every pair of clusters within each batch, analogous to the method used by [SingleR](https://www.bioconductor.org/packages/release/bioc/html/SingleR.html). In this case, we use the top 10 marker genes but any value can be used depending on the acceptable trade-off between signal and noise (and speed). We then take the union across all comparisons in all batches and use that in place of our HVG set in `fastMNN()`.


```r
# Recall that groups for marker detection
# are automatically defined from 'colLabels()'. 
markerList <- lapply(rescaled, function(x){
  y <- pairwiseWilcox(x, direction="up")
  getTopMarkers(y[[1]], y[[2]], n=10) %>% unlist %>% unlist
  })
marker.set <- unique(unlist(markerList))
#length(marker.set) # getting the total number of genes selected in this manner.
```

The total number of genes selected in this manner is: 426.


```r
set.seed(1000110)
mnn.out2 <- fastMNN(rescaled,
                    subset.row=marker.set,
                    BSPARAM=BiocSingular::RandomParam(deferred=TRUE))
# compute t-SNE:
mnn.out2 <- runTSNE(mnn.out2, dimred="corrected")
```

We can also visualize the corrected coordinates using a t-SNE plot:


```r
plotTSNE(mnn.out2, colour_by="batch")
```

<img src="dataSetIntegration_PRE-T_files/figure-html/unnamed-chunk-69-1.png" width="672" />

```r
plotTSNE(mnn.out2, colour_by="batch") + facet_wrap(~colData(mnn.out2)$batch)
```

<img src="dataSetIntegration_PRE-T_files/figure-html/unnamed-chunk-69-2.png" width="672" />

A quick inspection indicates that the original within-batch structure is indeed preserved in the corrected data. This highlights the utility of a marker-based feature set for integrating datasets that have already been characterized separately in a manner that preserves existing interpretations of each dataset. We note that some within-batch clusters have merged, most likely due to the lack of robust separation in the first place, though this may also be treated as a diagnostic on the appropriateness of the integration depending on the context.


```r
plotList <- vector(mode = "list", length = length(splVec))
for (x in 1:length(splVec)) {
  plotList[[x]] <- plotTSNE(mnn.out2[,mnn.out2$batch==splVec[x]],
                              colour_by=I(colLabels(rescaled[[x]]))) +
                  ggtitle(splVec[x])
}
gridExtra::grid.arrange(grobs = plotList,
      ncol=2
)
```

<img src="dataSetIntegration_PRE-T_files/figure-html/unnamed-chunk-70-1.png" width="672" />


```r
# by_exprs_values, : cannot find 'ENSG00000090382' for ETV6-RUNX1
# B cells
genex <- ensToShow[1]
p <- plotTSNE(mnn.out2, colour_by = genex, by_exprs_values="reconstructed")
p <- p + ggtitle(
  		paste("B cells", genex,
			rowData(uncorrected)[genex,"Symbol"])
			)
print(p)
```

```r
ensToShowPostCor <- ensToShow[ensToShow %in% rownames(mnn.out2)]
rowData(uncorrected)[ensToShowPostCor,c("ensembl_gene_id", "external_gene_name")]
```


```r
genex <- "ENSG00000156738" # MS4A1 # ensToShowPostCor[1]
p <- plotTSNE(mnn.out2, colour_by = genex, by_exprs_values="reconstructed")
p <- p + ggtitle(
  		paste(genex,
  		      rowData(uncorrected)[genex,"Symbol"])
			)
print(p)
```


```r
# ENSG00000203747 is FCGR3A
genex <- "ENSG00000203747" # FCGR3A # ensToShowPostCor[1]
p <- plotTSNE(mnn.out2, colour_by = genex, by_exprs_values="reconstructed")
p <- p + ggtitle(
  		paste(genex,
  		      rowData(uncorrected)[genex,"Symbol"])
			)
print(p)
```



```r
m.out <- findMarkers(uncorrected,
                     clusters.mnn,
                     block=uncorrected$batch,
                     direction="up",
                     lfc=1,
                     row.data=rowData(uncorrected)[,c("ensembl_gene_id","Symbol"),drop=FALSE])

#lapply(m.out, function(x){head(x[,2:6])})

# A (probably activated?) T cell subtype of some sort:
tl1 <- lapply(m.out, function(x){x[x$Symbol=="CD3D" & x$Top <= 50 & x$FDR < 0.10,2:6]}) # T-cell
tl2 <- lapply(m.out, function(x){x[x$Symbol=="CD69" & x$Top <= 50 & x$FDR < 0.20,2:6]}) # activation

tb1 <- unlist(lapply(tl1, nrow)) > 0
tb2 <- unlist(lapply(tl2, nrow)) > 0

cluToGet <- unique(c(which(tb1), which(tb2)))[1] # 3 # 19 # 4
demo <- m.out[[cluToGet]]
#as.data.frame(demo[1:20,c("Symbol", "Top", "p.value", "FDR", "summary.logFC")]) 
```

Expression level for the top gene,  on violin plots:


```r
geneEnsId <- rownames(demo)[1]
plotExpression(uncorrected, x=I(factor(clusters.mnn)), 
    features=geneEnsId, colour_by="batch") + facet_wrap(~colour_by) +
  ggtitle(sprintf("%s %s",
          geneEnsId,
          rowData(uncorrected)[geneEnsId,"Symbol"])
          )
```

<img src="dataSetIntegration_PRE-T_files/figure-html/unnamed-chunk-76-1.png" width="672" />

Expression level for the top gene, ENSG00000117632 on t-SNE plot:


```r
genex <- rownames(demo)[1]
p <- plotTSNE(mnn.out, colour_by = genex, by_exprs_values="reconstructed")
p <- p + ggtitle(
			paste("cluster", cluToGet, genex,
			rowData(uncorrected)[genex,"Symbol"])
		)
print(p)
```

<img src="dataSetIntegration_PRE-T_files/figure-html/unnamed-chunk-77-1.png" width="672" />


```r
p + facet_wrap(~colData(mnn.out)$batch)
```

<img src="dataSetIntegration_PRE-T_files/figure-html/unnamed-chunk-78-1.png" width="672" />

We suggest limiting the use of per-gene corrected values to visualization, e.g., when coloring points on a t-SNE plot by per-cell expression. This can be more aesthetically pleasing than uncorrected expression values that may contain large shifts on the colour scale between cells in different batches. Use of the corrected values in any quantitative procedure should be treated with caution, and should be backed up by similar results from an analysis on the uncorrected values.


```r
# save object?
#fn <- sprintf("dataSetIntegration_%s.Rdata", splSetToGet)
fn <- sprintf("%s/%s/Robjects/%s_sce_nz_postDeconv%s_dsi_%s.Rdata",
              projDir,
              outDirBit,
              setName,
              setSuf,
              splSetToGet) # 'dsi' for data set integration
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
##  [1] BiocSingular_1.6.0          Cairo_1.5-12.2             
##  [3] clustree_0.4.3              ggraph_2.0.5               
##  [5] pheatmap_1.0.12             dplyr_1.0.5                
##  [7] bluster_1.0.0               batchelor_1.6.3            
##  [9] scran_1.18.7                scater_1.18.6              
## [11] SingleCellExperiment_1.12.0 SummarizedExperiment_1.20.0
## [13] Biobase_2.50.0              GenomicRanges_1.42.0       
## [15] GenomeInfoDb_1.26.7         IRanges_2.24.1             
## [17] S4Vectors_0.28.1            BiocGenerics_0.36.1        
## [19] MatrixGenerics_1.2.1        matrixStats_0.58.0         
## [21] ggplot2_3.3.3               knitr_1.32                 
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-7              RColorBrewer_1.1-2       
##  [3] backports_1.2.1           tools_4.0.3              
##  [5] bslib_0.2.4               utf8_1.2.1               
##  [7] R6_2.5.0                  irlba_2.3.3              
##  [9] ResidualMatrix_1.0.0      vipor_0.4.5              
## [11] uwot_0.1.10               DBI_1.1.1                
## [13] colorspace_2.0-0          withr_2.4.2              
## [15] tidyselect_1.1.1          gridExtra_2.3            
## [17] compiler_4.0.3            BiocNeighbors_1.8.2      
## [19] DelayedArray_0.16.3       labeling_0.4.2           
## [21] bookdown_0.22             sass_0.3.1               
## [23] checkmate_2.0.0           scales_1.1.1             
## [25] stringr_1.4.0             digest_0.6.27            
## [27] rmarkdown_2.7             XVector_0.30.0           
## [29] pkgconfig_2.0.3           htmltools_0.5.1.1        
## [31] sparseMatrixStats_1.2.1   highr_0.9                
## [33] limma_3.46.0              rlang_0.4.10             
## [35] FNN_1.1.3                 DelayedMatrixStats_1.12.3
## [37] farver_2.1.0              jquerylib_0.1.3          
## [39] generics_0.1.0            jsonlite_1.7.2           
## [41] BiocParallel_1.24.1       RCurl_1.98-1.3           
## [43] magrittr_2.0.1            GenomeInfoDbData_1.2.4   
## [45] scuttle_1.0.4             Matrix_1.3-2             
## [47] Rcpp_1.0.6                ggbeeswarm_0.6.0         
## [49] munsell_0.5.0             fansi_0.4.2              
## [51] viridis_0.6.0             lifecycle_1.0.0          
## [53] stringi_1.5.3             yaml_2.2.1               
## [55] edgeR_3.32.1              MASS_7.3-54              
## [57] zlibbioc_1.36.0           Rtsne_0.15               
## [59] grid_4.0.3                ggrepel_0.9.1            
## [61] dqrng_0.3.0               crayon_1.4.1             
## [63] lattice_0.20-44           cowplot_1.1.1            
## [65] graphlayouts_0.7.1        beachmat_2.6.4           
## [67] locfit_1.5-9.4            pillar_1.6.0             
## [69] igraph_1.2.6              codetools_0.2-18         
## [71] glue_1.4.2                evaluate_0.14            
## [73] tweenr_1.0.2              vctrs_0.3.7              
## [75] polyclip_1.10-0           gtable_0.3.0             
## [77] purrr_0.3.4               tidyr_1.1.3              
## [79] assertthat_0.2.1          ggforce_0.3.3            
## [81] xfun_0.22                 rsvd_1.0.5               
## [83] tidygraph_1.2.0           RSpectra_0.16-0          
## [85] viridisLite_0.4.0         tibble_3.1.1             
## [87] beeswarm_0.3.1            statmod_1.4.35           
## [89] ellipsis_0.3.2
```

