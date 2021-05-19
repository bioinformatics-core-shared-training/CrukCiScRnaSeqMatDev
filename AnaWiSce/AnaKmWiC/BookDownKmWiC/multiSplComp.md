---
title: "CRUK CI Summer School 2020 - introduction to single-cell RNA-seq analysis"
subtitle: 'Multi-sample comparisons'

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
  splSetToGet: "PBMMC,ETV6-RUNX1"
  setName: "caron"
  setSuf: "_allCells"
  dsiSuf: '_dsi'
---

<!--
  setSuf: "_5hCellPerSpl"
  setSuf: "_5hCps"
-->
  

```r
projDir <- params$projDir
dirRel <- params$dirRel
outDirBit <- params$outDirBit
cacheBool <- params$cacheBool
splSetToGet <- params$splSetToGet
setName <- params$setName
setSuf <- params$setSuf
dsiSuf <- params$dsiSuf # 'dsi' for data set integration

if(params$bookType == "mk"){
	setName <- "caron"
	splSetToGet <- "PBMMC,ETV6-RUNX1"
	setSuf <- "_5hCps"
}

splSetVec <- unlist(strsplit(splSetToGet, ",")) # params may not be read in if knitting book.
splSetToGet2 <- gsub(",", "_", splSetToGet)
nbPcToComp <- 50
figSize <- 7
```





# Differential expression and abundance between conditions

Source: [Multi-sample comparisons](https://osca.bioconductor.org/multi-sample-comparisons.html) of the OSCA book.

## Motivation

A powerful use of scRNA-seq technology lies in the design of replicated multi-condition experiments to detect changes in composition or expression between conditions. For example, a researcher could use this strategy to detect changes in cell type abundance after drug treatment (Richard et al. 2018) or genetic modifications (Scialdone et al. 2016). This provides more biological insight than conventional scRNA-seq experiments involving only one biological condition, especially if we can relate population changes to specific experimental perturbations.

Differential analyses of multi-condition scRNA-seq experiments can be broadly split into two categories - differential expression (DE) and differential abundance (DA) analyses. The former tests for changes in expression between conditions for cells of the same type that are present in both conditions, while the latter tests for changes in the composition of cell types (or states, etc.) between conditions.

## Setting up the data

We will use the data set comprising the 11 samples (500 or 1000 cells per sample) analysed with fastMNN and the nested list of samples.

The differential analyses in this chapter will be predicated on many of the pre-processing steps covered previously. For brevity, we will not explicitly repeat them here, only noting that we have already merged cells from all samples into the same coordinate system and clustered the merged dataset to obtain a common partitioning across all samples.

Load the SCE object:


```r
#setName <- "caron"
# Read object in:
##setSuf <- "_1kCellPerSpl"
##tmpFn <- sprintf("%s/%s/Robjects/%s_sce_nz_postDeconv%s_clustered.Rds", projDir, outDirBit, setName, setSuf)

#setSuf <- "_1kCps"
tmpFn <- sprintf("%s/%s/Robjects/%s_sce_nz_postDeconv%s_Fmwbl.Rds", projDir, outDirBit, setName, setSuf)

print(tmpFn)
```

```
## [1] "/ssd/personal/baller01/20200511_FernandesM_ME_crukBiSs2020/AnaWiSce/AnaKmWiC/Robjects/caron_sce_nz_postDeconv_allCells_Fmwbl.Rds"
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
## dim: 12466 47830 
## metadata(2): merge.info pca.info
## assays(1): reconstructed
## rownames(12466): ENSG00000000003 ENSG00000000457 ... ENSG00000285476
##   ENSG00000285492
## rowData names(1): rotation
## colnames: NULL
## colData names(22): Sample Barcode ... type clusters.mnn
## reducedDimNames(2): corrected TSNE
## altExpNames(0):
```

A brief inspection of the results shows clusters contain varying contributions from batches:


```r
library(scater)
colLabels(sce) <- sce$clusters.mnn
tab <- table(colLabels(sce), sce$type)
tab
```

```
##      
##       ETV6-RUNX1  HHD PBMMC PRE-T
##   c1           3    2     6   904
##   c10       2110 1019   102    21
##   c11        212    7   235    17
##   c12         88    1     6     0
##   c13        385   92   112    17
##   c14         89    2    20     1
##   c15         37    1   143    21
##   c16         59    4   138    14
##   c17        164    3   125     6
##   c18       1745  295  2169   372
##   c19        649  418    65    30
##   c2        5272 2952   190    12
##   c20        125   81   410   548
##   c21          6    1    55    17
##   c22         30   18   318    49
##   c23         38    0     7     0
##   c24        475  206    42    34
##   c25        231  466   141    22
##   c26          1    2    23  1372
##   c27          2   16    85   348
##   c28          5    4   491    89
##   c29         16    2    48    12
##   c3          18   58    76   760
##   c30          3    1   270     2
##   c31         41    0    70     0
##   c32          4    0   199    25
##   c33         59    4   140     0
##   c34          0    0    80    15
##   c35          4    1    20     3
##   c36       4250 2137   203     6
##   c37        164    4   391     9
##   c38         83   52   100    28
##   c39        363  167   529    75
##   c4         798  575    63   998
##   c40         19    6    33    14
##   c41          0    1     7   554
##   c42        532   43     3     0
##   c43         10    2    26     4
##   c44          1   15    30   260
##   c45        329   96   327    62
##   c46          1    2    74     5
##   c47          3    1   140    52
##   c48          5    3    18     3
##   c49          3    4    20   100
##   c5          13    0   127     9
##   c50        117    9   171    15
##   c51        161    1   217     1
##   c52          1    2    99     3
##   c53         27    1    24     4
##   c54         49    0     0     0
##   c55          2    0    48     1
##   c56         20    0    26     0
##   c57         55   54   103    33
##   c58          7   19    78     6
##   c59          8   16    45   184
##   c6         414  167   457   223
##   c60          2    1    23     3
##   c61         40    0     1     0
##   c62          2    0   439     2
##   c63          1    1    80     4
##   c64          0    0    26     0
##   c65         26    0     0     1
##   c66          8    1    30     0
##   c7          40   51   352    67
##   c8          61  900   579     4
##   c9           2    1   210    28
```



```r
pheatmap::pheatmap(tab,
           border_color      = NA,
           drop_levels       = TRUE,
           cluster_rows      = FALSE,
           cluster_cols      = FALSE
           )
```

<img src="multiSplComp_files/figure-html/unnamed-chunk-5-1.png" width="672" />


```r
tab <- table(colLabels(sce), sce$Sample.Name2)
pheatmap::pheatmap(tab,
           border_color      = NA,
           drop_levels       = TRUE,
           cluster_rows      = FALSE,
           cluster_cols      = FALSE
           )
```

<img src="multiSplComp_files/figure-html/unnamed-chunk-6-1.png" width="672" />

On the t-SNE plots below, cells colored by type or sample ('batch of origin'). Cluster numbers are superimposed based on the median coordinate of cells assigned to that cluster. 


```r
p1 <- plotTSNE(sce, colour_by="type", text_by="label")
p2 <- plotTSNE(sce, colour_by="Sample.Name2")
gridExtra::grid.arrange(p1, p2+facet_wrap(~colData(sce)$type), ncol=2)
```

<img src="multiSplComp_files/figure-html/unnamed-chunk-7-1.png" width="806.4" />



```r
tmpFn <- sprintf("%s/%s/Robjects/%s_sce_nz_postDeconv%s_Fmwbl2.Rds", projDir, outDirBit, setName, setSuf)
tmpList <- readRDS(tmpFn)

chosen.hvgs <- tmpList$chosen.hvgs
rescaled.mbn <- tmpList$rescaled.mbn
uncorrected <- tmpList$uncorrected
colToKeep <- c("Run", "Sample.Name", "source_name", "block", "setName", "Sample.Name2") 
colData(uncorrected) <- colData(uncorrected)[,colToKeep]
colData(uncorrected)[1:3,]
```

```
## DataFrame with 3 rows and 6 columns
##           Run Sample.Name source_name      block     setName Sample.Name2
##   <character> <character>    <factor>   <factor> <character>  <character>
## 1  SRR9264343  GSM3872434  ETV6-RUNX1 ETV6-RUNX1       Caron ETV6-RUNX1_1
## 2  SRR9264343  GSM3872434  ETV6-RUNX1 ETV6-RUNX1       Caron ETV6-RUNX1_1
## 3  SRR9264343  GSM3872434  ETV6-RUNX1 ETV6-RUNX1       Caron ETV6-RUNX1_1
```

```r
#--- merging ---#
library(batchelor)
set.seed(01001001)
merged <- correctExperiments(uncorrected, 
    batch=uncorrected$Sample.Name2, 
    subset.row=chosen.hvgs,
    PARAM=FastMnnParam(
        merge.order=list( list(1,2,3,4), list(9,10,11), list(5,6), list(7,8) )
    )
)

merged
```

```
## class: SingleCellExperiment 
## dim: 12466 47830 
## metadata(2): merge.info pca.info
## assays(3): reconstructed counts logcounts
## rownames(12466): ENSG00000000003 ENSG00000000457 ... ENSG00000285476
##   ENSG00000285492
## rowData names(12): rotation ensembl_gene_id ... detected gene_sparsity
## colnames: NULL
## colData names(7): batch Run ... setName Sample.Name2
## reducedDimNames(4): corrected PCA TSNE UMAP
## altExpNames(0):
```

```r
#--- clustering ---#
g <- buildSNNGraph(merged, use.dimred="corrected")
clusters <- igraph::cluster_louvain(g)
merged$clusters.mnn <- factor(paste0("c", clusters$membership))
#colLabels(merged) <- merged$clusters.mnn

#--- dimensionality-reduction ---#
merged <- runTSNE(merged, dimred="corrected", external_neighbors=TRUE)
merged <- runUMAP(merged, dimred="corrected", external_neighbors=TRUE)

library(scater)
tab <- table(merged$clusters.mnn, merged$block)
pheatmap::pheatmap(tab,
           border_color      = NA,
           drop_levels       = TRUE,
           cluster_rows      = FALSE,
           cluster_cols      = FALSE
           )
```

<img src="multiSplComp_files/figure-html/unnamed-chunk-8-1.png" width="672" />

```r
tab <- table(merged$clusters.mnn, merged$Sample.Name2)
pheatmap::pheatmap(tab,
           border_color      = NA,
           drop_levels       = TRUE,
           cluster_rows      = FALSE,
           cluster_cols      = FALSE
           )
```

<img src="multiSplComp_files/figure-html/unnamed-chunk-8-2.png" width="672" />

```r
#plotTSNE(merged, colour_by="block", text_by="clusters.mnn")
#plotTSNE(merged, colour_by="Sample.Name2")
```


```r
p1 <- plotTSNE(merged, colour_by="block", text_by="clusters.mnn")
p2 <- plotTSNE(merged, colour_by="Sample.Name2")
gridExtra::grid.arrange(p1, p2+facet_wrap(~colData(sce)$type), ncol=2)
```

<img src="multiSplComp_files/figure-html/unnamed-chunk-9-1.png" width="806.4" />

## Differential expression between conditions

### Creating pseudo-bulk samples

The most obvious differential analysis is to look for changes in expression between conditions. We perform the DE analysis separately for each label. The actual DE testing is performed on “pseudo-bulk” expression profiles (Tung et al. 2017), generated by summing counts together for all cells with the same combination of label and sample. This leverages the resolution offered by single-cell technologies to define the labels, and combines it with the statistical rigor of existing methods for DE analyses involving a small number of samples.



```r
# Using 'label' and 'sample' as our two factors; each column of the output
# corresponds to one unique combination of these two factors.
summed <- aggregateAcrossCells(merged, 
    				id = DataFrame(
    					label=merged$clusters.mnn,
    					sample=merged$Sample.Name2
					)
)
summed
```

```
## class: SingleCellExperiment 
## dim: 12466 187 
## metadata(2): merge.info pca.info
## assays(1): counts
## rownames(12466): ENSG00000000003 ENSG00000000457 ... ENSG00000285476
##   ENSG00000285492
## rowData names(12): rotation ensembl_gene_id ... detected gene_sparsity
## colnames: NULL
## colData names(11): batch Run ... sample ncells
## reducedDimNames(4): corrected PCA TSNE UMAP
## altExpNames(0):
```

```r
colData(summed) %>% head(3)
```

```
## DataFrame with 3 rows and 11 columns
##          batch         Run Sample.Name source_name      block     setName
##    <character> <character> <character>    <factor>   <factor> <character>
## 1 ETV6-RUNX1_1  SRR9264343  GSM3872434  ETV6-RUNX1 ETV6-RUNX1       Caron
## 2 ETV6-RUNX1_2  SRR9264344  GSM3872435  ETV6-RUNX1 ETV6-RUNX1       Caron
## 3 ETV6-RUNX1_3  SRR9264345  GSM3872436  ETV6-RUNX1 ETV6-RUNX1       Caron
##   Sample.Name2 clusters.mnn    label       sample    ncells
##    <character>     <factor> <factor>  <character> <integer>
## 1 ETV6-RUNX1_1           c1       c1 ETV6-RUNX1_1       937
## 2 ETV6-RUNX1_2           c1       c1 ETV6-RUNX1_2       722
## 3 ETV6-RUNX1_3           c1       c1 ETV6-RUNX1_3       165
```

At this point, it is worth reflecting on the motivations behind the use of pseudo-bulking:

Larger counts are more amenable to standard DE analysis pipelines designed for bulk RNA-seq data. Normalization is more straightforward and certain statistical approximations are more accurate e.g., the saddlepoint approximation for quasi-likelihood methods or normality for linear models.
Collapsing cells into samples reflects the fact that our biological replication occurs at the sample level (Lun and Marioni 2017). Each sample is represented no more than once for each condition, avoiding problems from unmodelled correlations between samples. Supplying the per-cell counts directly to a DE analysis pipeline would imply that each cell is an independent biological replicate, which is not true from an experimental perspective. (A mixed effects model can handle this variance structure but involves extra statistical and computational complexity for little benefit, see Crowell et al. (2019).)
Variance between cells within each sample is masked, provided it does not affect variance across (replicate) samples. This avoids penalizing DEGs that are not uniformly up- or down-regulated for all cells in all samples of one condition. Masking is generally desirable as DEGs - unlike marker genes - do not need to have low within-sample variance to be interesting, e.g., if the treatment effect is consistent across replicate populations but heterogeneous on a per-cell basis. (Of course, high per-cell variability will still result in weaker DE if it affects the variability across populations, while homogeneous per-cell responses will result in stronger DE due to a larger population-level log-fold change. These effects are also largely desirable.)



### Performing the DE analysis

#### Introduction

The DE analysis will be performed using quasi-likelihood (QL) methods from the edgeR package (Robinson, McCarthy, and Smyth 2010; Chen, Lun, and Smyth 2016). This uses a negative binomial generalized linear model (NB GLM) to handle overdispersed count data in experiments with limited replication. In our case, we have biological variation with three paired replicates per condition, so edgeR (or its contemporaries) is a natural choice for the analysis.

We do not use all labels for GLM fitting as the strong DE between labels makes it difficult to compute a sensible average abundance to model the mean-dispersion trend. Moreover, label-specific batch effects would not be easily handled with a single additive term in the design matrix for the batch. Instead, we arbitrarily pick one of the labels to use for this demonstration.


```r
labelToGet <- "c1"
current <- summed[,summed$label==labelToGet]

# Creating up a DGEList object for use in edgeR:
suppressMessages(library(edgeR))
y <- DGEList(counts(current), samples=colData(current))
y
```

```
## An object of class "DGEList"
## $counts
##                 Sample1 Sample2 Sample3 Sample4 Sample5 Sample6 Sample7 Sample8
## ENSG00000000003       0       1       0       0       0       2       0       0
## ENSG00000000457      25      24       2      43       4       6       0       1
## ENSG00000000938       0       0       0       0       2       0       0       1
## ENSG00000001167      31      43      13      87       1      29       2       3
## ENSG00000001461      71      16       2      37      15      31       0       1
##                 Sample9 Sample10 Sample11
## ENSG00000000003       0        2        0
## ENSG00000000457       0        0        0
## ENSG00000000938       1        1        0
## ENSG00000001167       4        5        0
## ENSG00000001461       0        1        0
## 12461 more rows ...
## 
## $samples
##         group lib.size norm.factors        batch        Run Sample.Name
## Sample1     1  2399260            1 ETV6-RUNX1_1 SRR9264343  GSM3872434
## Sample2     1  1158608            1 ETV6-RUNX1_2 SRR9264344  GSM3872435
## Sample3     1   289586            1 ETV6-RUNX1_3 SRR9264345  GSM3872436
## Sample4     1  2161489            1 ETV6-RUNX1_4 SRR9264346  GSM3872437
## Sample5     1   773606            1        HHD_1 SRR9264347  GSM3872438
##         source_name      block setName Sample.Name2 clusters.mnn label
## Sample1  ETV6-RUNX1 ETV6-RUNX1   Caron ETV6-RUNX1_1           c1    c1
## Sample2  ETV6-RUNX1 ETV6-RUNX1   Caron ETV6-RUNX1_2           c1    c1
## Sample3  ETV6-RUNX1 ETV6-RUNX1   Caron ETV6-RUNX1_3           c1    c1
## Sample4  ETV6-RUNX1 ETV6-RUNX1   Caron ETV6-RUNX1_4           c1    c1
## Sample5         HHD        HHD   Caron        HHD_1           c1    c1
##               sample ncells
## Sample1 ETV6-RUNX1_1    937
## Sample2 ETV6-RUNX1_2    722
## Sample3 ETV6-RUNX1_3    165
## Sample4 ETV6-RUNX1_4   1363
## Sample5        HHD_1    276
## 6 more rows ...
```

#### Pre-processing

A typical step in bulk RNA-seq data analyses is to remove samples with very low library sizes due to failed library preparation or sequencing. The very low counts in these samples can be troublesome in downstream steps such as normalization (Chapter 7) or for some statistical approximations used in the DE analysis. In our situation, this is equivalent to removing label-sample combinations that have very few or lowly-sequenced cells. The exact definition of “very low” will vary, but in this case, we remove combinations containing fewer than 20 cells (Crowell et al. 2019). Alternatively, we could apply the outlier-based strategy described in Chapter 6, but this makes the strong assumption that all label-sample combinations have similar numbers of cells that are sequenced to similar depth.

<!--
with 500 cells per samples, some clusters are discarded, which may not be with more cells
-->


```r
discarded <- current$ncells < 20
y <- y[,!discarded]
summary(discarded)
```

```
##    Mode   FALSE    TRUE 
## logical      10       1
```

Another typical step in bulk RNA-seq analyses is to remove genes that are lowly expressed. This reduces computational work, improves the accuracy of mean-variance trend modelling and decreases the severity of the multiple testing correction. Genes are discarded if they are not expressed above a log-CPM threshold in a minimum number of samples (determined from the size of the smallest treatment group in the experimental design).


```r
keep <- filterByExpr(y, group=current$source_name)
y <- y[keep,]
summary(keep)
```

```
##    Mode   FALSE    TRUE 
## logical    6478    5988
```

Finally, we correct for composition biases by computing normalization factors with the trimmed mean of M-values method (Robinson and Oshlack 2010). We do not need the bespoke single-cell methods described in Chapter 7, as the counts for our pseudo-bulk samples are large enough to apply bulk normalization methods. (Readers should be aware that edgeR normalization factors are closely related but not the same as the size factors described elsewhere in this book.)


```r
y <- calcNormFactors(y)
y$samples
```

```
##          group lib.size norm.factors        batch        Run Sample.Name
## Sample1      1  2399260    0.7743944 ETV6-RUNX1_1 SRR9264343  GSM3872434
## Sample2      1  1158608    0.9666536 ETV6-RUNX1_2 SRR9264344  GSM3872435
## Sample3      1   289586    0.9940259 ETV6-RUNX1_3 SRR9264345  GSM3872436
## Sample4      1  2161489    0.8767754 ETV6-RUNX1_4 SRR9264346  GSM3872437
## Sample5      1   773606    0.6914457        HHD_1 SRR9264347  GSM3872438
## Sample6      1  1117032    0.8861710        HHD_2 SRR9264348  GSM3872439
## Sample7      1    37637    1.1549506      PBMMC_1 SRR9264351  GSM3872442
## Sample8      1    45151    1.3305647      PBMMC_2 SRR9264353  GSM3872443
## Sample9      1    37099    1.2773938      PBMMC_3 SRR9264354  GSM3872444
## Sample10     1    66418    1.2743275      PRE-T_1 SRR9264349  GSM3872440
##          source_name      block setName Sample.Name2 clusters.mnn label
## Sample1   ETV6-RUNX1 ETV6-RUNX1   Caron ETV6-RUNX1_1           c1    c1
## Sample2   ETV6-RUNX1 ETV6-RUNX1   Caron ETV6-RUNX1_2           c1    c1
## Sample3   ETV6-RUNX1 ETV6-RUNX1   Caron ETV6-RUNX1_3           c1    c1
## Sample4   ETV6-RUNX1 ETV6-RUNX1   Caron ETV6-RUNX1_4           c1    c1
## Sample5          HHD        HHD   Caron        HHD_1           c1    c1
## Sample6          HHD        HHD   Caron        HHD_2           c1    c1
## Sample7        PBMMC      PBMMC   Caron      PBMMC_1           c1    c1
## Sample8        PBMMC      PBMMC   Caron      PBMMC_2           c1    c1
## Sample9        PBMMC      PBMMC   Caron      PBMMC_3           c1    c1
## Sample10       PRE-T      PRE-T   Caron      PRE-T_1           c1    c1
##                sample ncells
## Sample1  ETV6-RUNX1_1    937
## Sample2  ETV6-RUNX1_2    722
## Sample3  ETV6-RUNX1_3    165
## Sample4  ETV6-RUNX1_4   1363
## Sample5         HHD_1    276
## Sample6         HHD_2    345
## Sample7       PBMMC_1     26
## Sample8       PBMMC_2     27
## Sample9       PBMMC_3     41
## Sample10      PRE-T_1     37
```

#### Statistical modelling

Our aim is to test whether the log-fold change between sample groups is significantly different from zero.


```r
design <- model.matrix(~factor(source_name), y$samples)
design
```

```
##          (Intercept) factor(source_name)HHD factor(source_name)PBMMC
## Sample1            1                      0                        0
## Sample2            1                      0                        0
## Sample3            1                      0                        0
## Sample4            1                      0                        0
## Sample5            1                      1                        0
## Sample6            1                      1                        0
## Sample7            1                      0                        1
## Sample8            1                      0                        1
## Sample9            1                      0                        1
## Sample10           1                      0                        0
##          factor(source_name)PRE-T
## Sample1                         0
## Sample2                         0
## Sample3                         0
## Sample4                         0
## Sample5                         0
## Sample6                         0
## Sample7                         0
## Sample8                         0
## Sample9                         0
## Sample10                        1
## attr(,"assign")
## [1] 0 1 1 1
## attr(,"contrasts")
## attr(,"contrasts")$`factor(source_name)`
## [1] "contr.treatment"
```

We estimate the negative binomial (NB) dispersions with estimateDisp(). The role of the NB dispersion is to model the mean-variance trend, which is not easily accommodated by QL dispersions alone due to the quadratic nature of the NB mean-variance trend.


```r
y <- estimateDisp(y, design)
summary(y$trended.dispersion)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.1498  0.1546  0.1639  0.1843  0.1876  0.3670
```

Biological coefficient of variation (BCV) for each gene as a function of the average abundance. The BCV is computed as the square root of the NB dispersion after empirical Bayes shrinkage towards the trend. Trended and common BCV estimates are shown in blue and red, respectively. 


```r
plotBCV(y)
```

<img src="multiSplComp_files/figure-html/unnamed-chunk-17-1.png" width="672" />

We also estimate the quasi-likelihood dispersions with glmQLFit() (Chen, Lun, and Smyth 2016). This fits a GLM to the counts for each gene and estimates the QL dispersion from the GLM deviance. We set robust=TRUE to avoid distortions from highly variable clusters (Phipson et al. 2016). The QL dispersion models the uncertainty and variability of the per-gene variance - which is not well handled by the NB dispersions, so the two dispersion types complement each other in the final analysis.


```r
fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$var.prior)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.4648  0.7190  0.8145  0.8089  0.8702  1.2078
```


```r
summary(fit$df.prior)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.3833 10.1200 10.1200  9.5004 10.1200 10.1200
```

QL dispersion estimates for each gene as a function of abundance. Raw estimates (black) are shrunk towards the trend (blue) to yield squeezed estimates (red).


```r
plotQLDisp(fit)
```

<img src="multiSplComp_files/figure-html/unnamed-chunk-20-1.png" width="672" />

We test for differences in expression due to sample group using glmQLFTest(). DEGs are defined as those with non-zero log-fold changes at a false discovery rate of 5%. If very few genes are significantly DE that sample group has little effect on the transcriptome.


```r
res <- glmQLFTest(fit, coef=ncol(design))
summary(decideTests(res))
```

```
##        factor(source_name)PRE-T
## Down                        186
## NotSig                     5726
## Up                           76
```


```r
topTab <- topTags(res)$table
tmpAnnot <- rowData(current)[,c("ensembl_gene_id","Symbol")] %>% data.frame
topTab %>% tibble::rownames_to_column("ensembl_gene_id") %>%
	left_join(tmpAnnot, by="ensembl_gene_id")
```

```
##    ensembl_gene_id      logFC    logCPM         F       PValue          FDR
## 1  ENSG00000137731  11.750497  8.081049 200.14071 1.137048e-09 6.808645e-06
## 2  ENSG00000100721 -14.574197 11.841139  86.06404 7.246895e-08 2.169720e-04
## 3  ENSG00000019582  -7.293701 13.419485  75.91130 1.699885e-07 3.392971e-04
## 4  ENSG00000030419   5.108326  4.324652  60.99556 7.199543e-07 9.130830e-04
## 5  ENSG00000204287  -5.991293 12.828441  60.45677 7.624274e-07 9.130830e-04
## 6  ENSG00000229989   3.928677  7.844280  57.24178 1.082802e-06 1.080637e-03
## 7  ENSG00000081189  -5.661665  9.823323  55.76542 1.278801e-06 1.093923e-03
## 8  ENSG00000272398  -7.136564 10.987138  53.63986 1.635052e-06 1.210049e-03
## 9  ENSG00000128218  -5.479618 11.067708  52.73912 1.818710e-06 1.210049e-03
## 10 ENSG00000100629   4.635353  4.820067  49.53548 2.687517e-06 1.422310e-03
##        Symbol
## 1       FXYD2
## 2       TCL1A
## 3        CD74
## 4       IKZF2
## 5     HLA-DRA
## 6  MIR181A1HG
## 7       MEF2C
## 8        CD24
## 9      VPREB3
## 10     CEP128
```

#### Differential expression for each cluster

The steps illustrated above with cluster c1 are now repeated for each cluster:

* Subset pseudo-bulk counts for that cluster
* Create edgeR object with these pseudo-bulk counts
* Pre-process
    * Remove samples with very small library size
    * Remove genes with low UMI counts
    * Correct for compositional bias
* Perform differential expression analysis  
    * Estimate negative binomial dispersion
    * Estimate quasi-likelihood dispersion
    * Test for differential expression 


```r
de.results <- list()
for (labelToGet in levels(summed$label)) {

	current <- summed[,summed$label==labelToGet]

    y <- DGEList(counts(current), samples=colData(current))

    discarded <- isOutlier(colSums(counts(current)), log=TRUE, type="lower")
    y <- y[,!discarded]
    y <- y[filterByExpr(y, group=current$source_name),]
    y <- calcNormFactors(y)

    design <- try(
        model.matrix(~factor(source_name), y$samples),
        silent=TRUE
    )
    if (is(design, "try-error") || 
        qr(design)$rank==nrow(design) ||
        qr(design)$rank < ncol(design)) 
    {
        # Skipping labels without contrasts or without 
        # enough residual d.f. to estimate the dispersion.
        next
    }

    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)
    res <- glmQLFTest(fit, coef=ncol(design))
    de.results[[labelToGet]] <- res
}
```

##### Number of DEGs by cluster and direction

We examine the numbers of DEGs at a FDR of 5% for each label (i.e. cluster). In general, there seems to be very little differential expression between the on and off conditions.


```r
summaries <- lapply(de.results, FUN=function(x) summary(decideTests(x))[,1])
sum.tab <- do.call(rbind, summaries)
#sum.tab
sum.tab[order(rownames(sum.tab)),] %>%
	as.data.frame() %>%
	tibble::rownames_to_column("Cluster") %>%
	datatable(rownames = FALSE, options = list(pageLength = 20, scrollX = TRUE))
```

```{=html}
<div id="htmlwidget-6e3d648fe7a73f057d76" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-6e3d648fe7a73f057d76">{"x":{"filter":"none","data":[["c1","c10","c11","c12","c13","c14","c15","c16","c17","c18","c2","c3","c4","c5","c6","c7","c8","c9"],[303,503,448,0,14,89,0,0,2,0,726,972,988,18,8,0,264,0],[3938,2532,2371,64,1255,2736,1868,544,1675,586,4010,5736,5352,4680,3685,229,1810,1380],[81,377,298,0,29,22,1,0,21,0,513,701,712,18,15,1,151,0]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>Cluster<\/th>\n      <th>Down<\/th>\n      <th>NotSig<\/th>\n      <th>Up<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"pageLength":20,"scrollX":true,"columnDefs":[{"className":"dt-right","targets":[1,2,3]}],"order":[],"autoWidth":false,"orderClasses":false,"lengthMenu":[10,20,25,50,100]}},"evals":[],"jsHooks":[]}</script>
```

##### List of DEGs

We now list DEGs and the number of clusters they were detected in:


```r
degs <- lapply(de.results, FUN=function(x) rownames(topTags(x, p.value=0.05)))
common.degs <- sort(table(unlist(degs)), decreasing=TRUE)
#head(common.degs, 20)
common.degs %>%
	as.data.frame %>% 
	dplyr::rename(ensembl_gene_id = Var1, NbClu = Freq) %>%
	left_join(
	data.frame(rowData(summed)[,c("ensembl_gene_id", "Symbol")]),
	by="ensembl_gene_id") %>%
	#rename(Gene = ensembl_gene_id) %>%
	relocate(c("Symbol","NbClu","ensembl_gene_id")) %>%
	datatable(rownames = FALSE, options = list(pageLength = 20, scrollX = TRUE))
```

```{=html}
<div id="htmlwidget-496041b9cf1cd0fed19d" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-496041b9cf1cd0fed19d">{"x":{"filter":"none","data":[["FXYD2","CHI3L2","MAL","MEF2C","TCL1A","S100A9","GNG7","HLA-DRA","VPREB3","S100A8","SH2D1A","CD24","CD74","CD79A","MDK","TNFAIP3","SOCS2","ALDH1A2","CD1E","SPON2","CD3G","CD3D","ARL4C","HLA-DMA","HLA-DPB1","HLA-DPA1","MPO","STAG3","TRAF4","TCF7","LYZ","NCF4","RETN","EPHB6","AHR","BIN2","KLRB1","CD2","CTGF","BCL11B","LILRB2","BEX1","RCBTB2","SLC27A3","HHEX","MYO7B","HOPX","AL138899.1","CD19","SH2D4B","KCTD12","ZNF703","PAX5","HLA-DRB1","PRTN3","ANKRD28","TRBC2","LAT","PPM1N","HLA-DQA2","AC007952.4","Z93241.1","AC103591.3"],[9,7,6,5,5,5,4,4,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],["ENSG00000137731","ENSG00000064886","ENSG00000172005","ENSG00000081189","ENSG00000100721","ENSG00000163220","ENSG00000176533","ENSG00000204287","ENSG00000128218","ENSG00000143546","ENSG00000183918","ENSG00000272398","ENSG00000019582","ENSG00000105369","ENSG00000110492","ENSG00000118503","ENSG00000120833","ENSG00000128918","ENSG00000158488","ENSG00000159674","ENSG00000160654","ENSG00000167286","ENSG00000188042","ENSG00000204257","ENSG00000223865","ENSG00000231389","ENSG00000005381","ENSG00000066923","ENSG00000076604","ENSG00000081059","ENSG00000090382","ENSG00000100365","ENSG00000104918","ENSG00000106123","ENSG00000106546","ENSG00000110934","ENSG00000111796","ENSG00000116824","ENSG00000118523","ENSG00000127152","ENSG00000131042","ENSG00000133169","ENSG00000136161","ENSG00000143554","ENSG00000152804","ENSG00000169994","ENSG00000171476","ENSG00000176320","ENSG00000177455","ENSG00000178217","ENSG00000178695","ENSG00000183779","ENSG00000196092","ENSG00000196126","ENSG00000196415","ENSG00000206560","ENSG00000211772","ENSG00000213658","ENSG00000213889","ENSG00000237541","ENSG00000262202","ENSG00000270022","ENSG00000273338"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>Symbol<\/th>\n      <th>NbClu<\/th>\n      <th>ensembl_gene_id<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"pageLength":20,"scrollX":true,"columnDefs":[{"className":"dt-right","targets":1}],"order":[],"autoWidth":false,"orderClasses":false,"lengthMenu":[10,20,25,50,100]}},"evals":[],"jsHooks":[]}</script>
```

##### Number of clusters skipped

"We also list the labels that were skipped due to the absence of replicates or contrasts. If it is necessary to extract statistics in the absence of replicates, several strategies can be applied such as reducing the complexity of the model or using a predefined value for the NB dispersion. We refer readers to the edgeR user’s guide for more details."


```r
skippedClusters <- setdiff(unique(summed$label), names(summaries))
```

The number of clusters skipped is 0.


```r
if(length(skippedClusters)>0)
{
  skippedClusters
}
```


```r
grmToShowList <- vector("list", length = nlevels(merged$clusters.mnn))
names(grmToShowList) <- levels(merged$clusters.mnn)
genesToExclude <- c()
nbGeneToShow <- 20

#degs <- lapply(de.results, FUN=function(x) (topTags(x, p.value=0.05)))
degs <- lapply(de.results, FUN=function(x) (as.data.frame(topTags(x, n=nbGeneToShow))))

for( namex in levels(merged$clusters.mnn) )
{
	nbGeneToUse <- min(c(nrow(degs[[namex]]), nbGeneToShow))

	# format

	# format p value:
	tmpCol <- grep("PValue|FDR", colnames(degs[[namex]]), value=TRUE)
	degs[[namex]][,tmpCol] <- apply(degs[[namex]][,tmpCol],
					     2,
					     function(x){format(x, scientific = TRUE, digits = 1)})
	# format logFC:
	tmpCol <- c("logFC", "logCPM", "F")
	degs[[namex]][,tmpCol] <- apply(degs[[namex]][,tmpCol], 2,  function(x){round(x, 2)})
	rm(tmpCol)

	# subset data
	grmToShow <- degs[[namex]] %>%
		as.data.frame() %>%
		tibble::rownames_to_column("gene") %>%	
		arrange(FDR, desc(abs(logFC))) %>%
		filter(! gene %in% genesToExclude) %>%
		group_modify(~ head(.x, nbGeneToUse)) 
	# keep data
	grmToShow$cluster <- namex
	grmToShowList[[namex]] <- grmToShow
	# tidy
	rm(nbGeneToUse)
}
grmToShowDf <- do.call("rbind", grmToShowList)
tmpCol <- c("cluster", "gene")
grmToShowDf %>%
	select(tmpCol, setdiff(colnames(grmToShowDf), tmpCol)) %>%
	filter(gene %in% names(common.degs) & as.numeric(FDR) < 0.05) %>%
	datatable(rownames = FALSE, filter="top", options=list(scrollX = TRUE, pageLength = 15))
```

```{=html}
<div id="htmlwidget-7177ff981b054aeb3d66" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-7177ff981b054aeb3d66">{"x":{"filter":"top","filterHTML":"<tr>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none; position: absolute; width: 200px;\">\n      <div data-min=\"-14.68\" data-max=\"14.3\" data-scale=\"2\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none; position: absolute; width: 200px;\">\n      <div data-min=\"3.1\" data-max=\"14.08\" data-scale=\"2\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none; position: absolute; width: 200px;\">\n      <div data-min=\"19\" data-max=\"475.64\" data-scale=\"2\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n<\/tr>","data":[["c1","c1","c1","c1","c1","c1","c1","c1","c1","c1","c1","c1","c1","c1","c1","c1","c1","c1","c1","c10","c10","c10","c10","c10","c10","c10","c10","c10","c10","c10","c10","c10","c10","c10","c11","c11","c11","c11","c11","c11","c11","c11","c11","c11","c11","c11","c11","c13","c13","c13","c13","c13","c13","c13","c13","c13","c13","c13","c14","c14","c14","c14","c14","c14","c14","c14","c14","c14","c14","c14","c14","c14","c14","c14","c14","c15","c17","c17","c17","c17","c17","c17","c17","c17","c17","c17","c2","c2","c2","c2","c2","c2","c2","c2","c2","c2","c2","c2","c2","c3","c3","c3","c3","c3","c3","c3","c3","c3","c3","c3","c3","c3","c4","c4","c4","c4","c4","c4","c4","c4","c4","c4","c4","c4","c4","c5","c5","c5","c5","c5","c5","c5","c5","c5","c5","c5","c5","c5","c5","c5","c5","c6","c6","c6","c6","c6","c6","c6","c6","c6","c6","c6","c6","c6","c6","c6","c6","c7","c8","c8","c8","c8","c8","c8","c8","c8","c8","c8","c8","c8","c8","c8","c8","c8","c8","c8"],["ENSG00000183918","ENSG00000272398","ENSG00000204287","ENSG00000213658","ENSG00000019582","ENSG00000176320","ENSG00000223865","ENSG00000116824","ENSG00000160654","ENSG00000172005","ENSG00000110492","ENSG00000231389","ENSG00000081189","ENSG00000128218","ENSG00000127152","ENSG00000100721","ENSG00000128918","ENSG00000064886","ENSG00000137731","ENSG00000231389","ENSG00000081189","ENSG00000137731","ENSG00000223865","ENSG00000116824","ENSG00000172005","ENSG00000159674","ENSG00000188042","ENSG00000110934","ENSG00000183918","ENSG00000064886","ENSG00000204257","ENSG00000213658","ENSG00000160654","ENSG00000127152","ENSG00000223865","ENSG00000100721","ENSG00000137731","ENSG00000118523","ENSG00000128218","ENSG00000152804","ENSG00000204287","ENSG00000081189","ENSG00000177455","ENSG00000064886","ENSG00000176533","ENSG00000183918","ENSG00000272398","ENSG00000196415","ENSG00000005381","ENSG00000167286","ENSG00000183779","ENSG00000100721","ENSG00000213889","ENSG00000104918","ENSG00000106546","ENSG00000118503","ENSG00000178695","ENSG00000163220","ENSG00000105369","ENSG00000231389","ENSG00000090382","ENSG00000137731","ENSG00000163220","ENSG00000143546","ENSG00000172005","ENSG00000100721","ENSG00000272398","ENSG00000158488","ENSG00000064886","ENSG00000128218","ENSG00000176533","ENSG00000204287","ENSG00000223865","ENSG00000273338","ENSG00000110492","ENSG00000163220","ENSG00000262202","ENSG00000143546","ENSG00000163220","ENSG00000237541","ENSG00000118503","ENSG00000270022","ENSG00000076604","ENSG00000273338","ENSG00000066923","ENSG00000171476","ENSG00000172005","ENSG00000064886","ENSG00000160654","ENSG00000128218","ENSG00000178217","ENSG00000196092","ENSG00000120833","ENSG00000159674","ENSG00000110492","ENSG00000131042","ENSG00000127152","ENSG00000143554","ENSG00000137731","ENSG00000158488","ENSG00000160654","ENSG00000064886","ENSG00000081189","ENSG00000137731","ENSG00000204287","ENSG00000172005","ENSG00000176533","ENSG00000081059","ENSG00000204257","ENSG00000133169","ENSG00000231389","ENSG00000100365","ENSG00000183918","ENSG00000167286","ENSG00000137731","ENSG00000160654","ENSG00000064886","ENSG00000128218","ENSG00000128918","ENSG00000176533","ENSG00000159674","ENSG00000110492","ENSG00000176320","ENSG00000133169","ENSG00000100721","ENSG00000272398","ENSG00000120833","ENSG00000105369","ENSG00000188042","ENSG00000106123","ENSG00000211772","ENSG00000100721","ENSG00000159674","ENSG00000128218","ENSG00000127152","ENSG00000137731","ENSG00000177455","ENSG00000171476","ENSG00000167286","ENSG00000064886","ENSG00000172005","ENSG00000196126","ENSG00000204287","ENSG00000206560","ENSG00000019582","ENSG00000137731","ENSG00000111796","ENSG00000106123","ENSG00000090382","ENSG00000143546","ENSG00000136161","ENSG00000110492","ENSG00000128218","ENSG00000272398","ENSG00000169994","ENSG00000100721","ENSG00000163220","ENSG00000163220","ENSG00000204257","ENSG00000204287","ENSG00000081189","ENSG00000137731","ENSG00000172005","ENSG00000019582","ENSG00000223865","ENSG00000105369","ENSG00000231389","ENSG00000176533","ENSG00000100365","ENSG00000064886","ENSG00000272398","ENSG00000128218","ENSG00000120833","ENSG00000196126","ENSG00000100721","ENSG00000128918"],[9.13,-7.35,-6.39,6.15,-7.36,10.62,-8.5,8.15,8.12,7.26,-7.06,-6.07,-5.87,-5.73,5.4,-14.45,11.71,9.93,11.14,-8.52,-12.68,10.73,-8.88,8.06,7.99,6.85,6.27,5.12,9.25,9.8,-8.52,6.94,9.25,6.92,-8.3,-10.4,11.26,-10.92,-10.09,-9.05,-8.75,-8.43,-8.25,7.92,-7.69,5.89,-8.24,11.3,10.6,9.57,-5.76,-4.69,-4.56,4.23,3.95,3.33,3.25,3.24,-5.94,-5.48,4.63,7.25,6.22,5.63,3.72,-10.47,-10.33,8.39,8.39,-6.69,-8.4,-5.83,-5.39,3.31,-6.52,7.53,4.11,6.1,6.01,4.08,3.17,2.87,2.36,4.09,3.25,4.17,14.3,10.68,10.53,-9.69,-9.67,-9.61,-8.65,7.11,-6.85,-6.77,6.56,-6.28,12.45,12.36,10.77,10.25,-10.24,10.22,-9.2,8.92,-8.48,5.25,-7.78,9.82,-9.18,-11.46,9.07,8.6,11.05,9.49,10.96,-11.81,11.26,-8.87,7.8,-7.39,9.64,9.36,-12.47,-13.83,-13.57,-6.05,5.33,5.18,4.95,-14.68,7.32,-6.27,4.88,9.81,-11.95,7.06,6.49,8.92,7.05,-4.48,-3.36,3.08,-2.17,4.69,3,2.1,7.86,5.88,3.2,-7.87,-4.24,-3.76,4.4,-8.85,5.04,5.93,-10.2,-8.9,-8.88,8.58,6.63,-7.36,-7.35,-6.41,-9.2,-8.58,-8.4,8.38,-6.85,-9.36,-8.37,-6.77,-8.15,7.21],[5.83,11.03,12.87,5.62,13.45,4.77,11.28,5.88,6.38,6.69,9.48,10.8,9.88,11.15,5.25,11.82,6.31,8.97,8.34,10.51,9.75,8.47,11.01,7.14,7.93,5.92,7.4,6.07,6.33,8.7,9.65,7.71,7.34,7.27,10.9,11.21,9.18,8.05,10.39,7.69,12.28,9.52,7.68,9.39,8.16,7.27,10.06,8.94,9.43,6.29,6.08,7.16,7.18,9.31,6.92,7.04,7.81,14.08,6.88,6.43,5.5,3.76,6.06,6.12,5.39,7.13,6.47,3.1,4.83,6.08,4.58,8.31,6.47,4.73,5.73,7.19,7.52,7.22,6.98,6.67,5.81,7.02,6.77,8.18,7.42,7.15,8.25,9.45,7.39,9.29,4.5,5.45,9.01,6.03,8.57,5.54,6.57,6.32,9.01,7.97,7.19,9.31,9.12,8.94,12.23,7.36,7.86,7.54,9.09,6.22,10.12,7.27,6.51,8.85,8.89,7.17,9.38,10.15,7.08,7.87,5.02,9.25,5.48,6.45,11.32,10.77,10.31,11.63,5.9,6.97,8.34,12.22,3.8,10.97,4.9,6.38,9.06,4.46,6.62,8.35,5.7,6.02,7.52,4.94,9.05,4.78,8.03,6.44,6.11,6.5,5.7,4.61,5.29,5.47,4.15,6.75,6.57,6.13,9.11,11.8,9.74,8.12,7.62,12.64,9.95,10.25,9.95,7.77,7.51,8.61,9.58,10.09,8.87,10.55,10.1,7.03],[150.85,67.43,68.65,74.47,106.15,105.86,59.95,61.91,57.53,69.46,57.13,56.89,60.75,60.11,62.46,98.09,200.65,114.6,145.07,135.27,124.22,109.33,101.49,110.46,98.11,113.13,101.77,107.94,202.68,128.87,91.37,92.03,76.55,76.86,83.16,145.74,111.19,111.93,106.33,104.78,109.2,109.6,110.55,103.59,103.12,110.75,94.53,25.44,25.93,26.89,28.23,24.73,26.44,23.46,22.56,23.13,25.06,20.12,28.66,30.27,29.17,50.95,48.4,46.38,47.86,42.41,43.18,40.92,76,35.86,35.03,33.94,33.07,32.04,31.27,37.88,40.79,25.21,24.9,27.49,22.53,25.82,23.32,58.69,52.3,47.37,239.11,166.32,124.93,130.26,156.85,207.1,140.06,130.51,142.45,139.16,140.07,144.15,241.43,331.69,297.12,225.1,212.18,245,239.82,214.08,221.45,217.71,195.68,202.9,164,188.26,146.89,152.43,318.92,275.52,475.64,216.31,300.2,193.99,215.15,199.75,269.38,184.34,165.2,36.44,37.49,33.96,36.06,37.57,36.34,29.64,32.45,31.63,31.08,167.77,25.72,27.7,72.06,107.77,97.49,25,24.12,25.66,24.21,62.92,35.56,20.85,55.64,47.83,48.5,19.31,19,19.46,30.31,29.44,41.53,31.38,95.02,110.21,87.18,102.34,92.74,73.86,76.03,73.29,69.55,70.97,75.88,78.35,66.3,60.08,60.89,59.82,57.56,65.98],["5e-09","3e-07","3e-07","5e-07","1e-08","9e-07","7e-07","1e-06","1e-06","1e-06","1e-06","1e-06","7e-07","7e-07","1e-06","5e-08","4e-08","3e-08","2e-09","1e-08","4e-08","5e-08","8e-08","4e-08","9e-08","4e-08","7e-08","5e-08","9e-10","2e-07","1e-07","1e-07","4e-07","4e-07","8e-08","1e-09","2e-08","2e-08","1e-08","2e-08","1e-08","1e-08","1e-08","2e-08","2e-08","1e-08","3e-08","2e-04","2e-04","1e-04","6e-05","1e-04","8e-05","2e-04","2e-04","2e-04","1e-04","3e-04","9e-06","6e-06","9e-06","1e-07","1e-07","2e-07","1e-07","4e-07","4e-07","9e-07","1e-09","1e-06","2e-06","2e-06","3e-06","4e-06","5e-06","1e-05","3e-06","7e-05","7e-05","5e-05","1e-04","6e-05","1e-04","2e-07","6e-07","1e-06","2e-08","8e-09","4e-08","3e-08","3e-08","5e-09","2e-08","3e-08","2e-08","2e-08","2e-08","2e-08","8e-10","8e-10","4e-10","2e-09","3e-09","1e-09","2e-09","3e-09","3e-09","3e-09","5e-09","1e-08","1e-08","2e-08","3e-08","2e-08","3e-10","8e-10","3e-11","3e-09","5e-09","6e-09","3e-09","5e-09","9e-09","8e-09","2e-08","2e-05","2e-05","3e-05","2e-05","3e-05","2e-05","7e-05","6e-05","4e-05","4e-05","4e-09","1e-04","1e-04","3e-07","2e-08","4e-08","4e-05","5e-05","3e-05","5e-05","3e-08","3e-06","1e-04","2e-07","3e-07","3e-07","2e-04","2e-04","2e-04","1e-05","1e-05","9e-07","6e-05","4e-08","8e-09","4e-08","4e-08","3e-08","1e-07","1e-07","1e-07","2e-07","3e-07","2e-07","3e-07","3e-07","6e-07","5e-07","6e-07","8e-07","7e-07"],["1e-05","2e-04","2e-04","2e-04","2e-05","3e-04","3e-04","3e-04","3e-04","3e-04","3e-04","3e-04","3e-04","3e-04","3e-04","3e-05","3e-05","3e-05","7e-06","2e-05","3e-05","3e-05","3e-05","3e-05","3e-05","3e-05","3e-05","3e-05","3e-06","4e-05","4e-05","4e-05","7e-05","7e-05","1e-05","4e-06","5e-06","5e-06","5e-06","5e-06","5e-06","5e-06","5e-06","5e-06","5e-06","5e-06","8e-06","2e-02","2e-02","2e-02","2e-02","2e-02","2e-02","2e-02","2e-02","2e-02","2e-02","3e-02","1e-03","1e-03","1e-03","1e-04","1e-04","1e-04","1e-04","2e-04","2e-04","3e-04","4e-06","5e-04","6e-04","6e-04","7e-04","8e-04","9e-04","2e-02","1e-03","2e-02","2e-02","2e-02","2e-02","2e-02","2e-02","4e-04","5e-04","6e-04","1e-05","1e-05","1e-05","1e-05","1e-05","1e-05","1e-05","1e-05","1e-05","1e-05","1e-05","1e-05","4e-06","3e-06","3e-06","3e-06","3e-06","3e-06","3e-06","3e-06","3e-06","3e-06","4e-06","6e-06","6e-06","7e-06","1e-05","1e-05","1e-06","2e-06","2e-07","5e-06","5e-06","5e-06","5e-06","5e-06","6e-06","6e-06","8e-06","1e-02","1e-02","1e-02","1e-02","1e-02","1e-02","2e-02","2e-02","2e-02","2e-02","2e-05","3e-02","3e-02","3e-04","4e-05","6e-05","1e-02","1e-02","1e-02","1e-02","1e-04","2e-03","3e-02","3e-04","3e-04","3e-04","4e-02","4e-02","4e-02","5e-03","7e-03","7e-04","1e-02","2e-05","2e-05","2e-05","2e-05","2e-05","4e-05","4e-05","4e-05","5e-05","5e-05","5e-05","5e-05","5e-05","7e-05","7e-05","7e-05","8e-05","8e-05"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>cluster<\/th>\n      <th>gene<\/th>\n      <th>logFC<\/th>\n      <th>logCPM<\/th>\n      <th>F<\/th>\n      <th>PValue<\/th>\n      <th>FDR<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"scrollX":true,"pageLength":15,"columnDefs":[{"className":"dt-right","targets":[2,3,4]}],"order":[],"autoWidth":false,"orderClasses":false,"orderCellsTop":true,"lengthMenu":[10,15,25,50,100]}},"evals":[],"jsHooks":[]}</script>
```

```r
tmpBool <- as.numeric(grmToShowDf$FDR) < 0.05 
markers.to.plot <- unique(grmToShowDf[tmpBool, "gene"])
markers.to.plot <- markers.to.plot[1:5]
```

### Putting it all together

Now that we have laid out the theory underlying the DE analysis, we repeat this process for each of the labels. This is conveniently done using the pseudoBulkDGE() function from scran, which will loop over all labels and apply the exact analysis described above to each label. To prepare for this, we filter out all sample-label combinations with insufficient cells.


```r
summed.filt <- summed[,summed$ncells >= 20]
```

We construct a common design matrix that will be used in the analysis for each label. Recall that this matrix should have one row per unique sample (and named as such), reflecting the fact that we are modelling counts on the sample level instead of the cell level.


```r
# Pulling out a sample-level 'targets' data.frame:
targets <- colData(merged)[!duplicated(merged$Sample.Name2),]

# Constructing the design matrix:
design <- model.matrix(~factor(source_name), data=targets)
rownames(design) <- targets$Sample.Name2
```

We then apply the pseudoBulkDGE() function to obtain a list of DE genes for each label. This function puts some additional effort into automatically dealing with labels that are not represented in all sample groups, for which a DE analysis between conditions is meaningless; or are not represented in a sufficient number of replicate samples to enable modelling of biological variability.


```r
library(scran)
de.results <- pseudoBulkDGE(summed.filt, 
    sample=summed.filt$Sample.Name2,
    label=summed.filt$label,
    design=design,
    coef=ncol(design),

    # 'condition' sets the group size for filterByExpr(),
    # to perfectly mimic our previous manual analysis.
    condition=targets$source_name 
)
```

We examine the numbers of DEGs at a FDR of 5% for each label using the decideTestsPerLabel() function. Note that genes listed as NA were either filtered out as low-abundance genes for a given label’s analysis, or the comparison of interest was not possible for a particular label, e.g., due to lack of residual degrees of freedom or an absence of samples from both conditions.


```r
is.de <- decideTestsPerLabel(de.results, threshold=0.05)
summarizeTestsPerLabel(is.de)
```

```
##      -1    0   1    NA
## c1  202 6786  91  5387
## c10 621 4186 429  7230
## c11 462 2763 355  8886
## c13  64 5309 206  6887
## c15   0 2829   0  9637
## c17   5 2056  20 10385
## c2  708 4033 508  7217
## c3  942 5766 701  5057
## c4  951 5396 705  5414
## c6    6 3701  12  8747
## c8  257 1764 156 10289
```

For each gene, we compute the percentage of cell types in which that gene is upregulated or downregulated. (Here, we consider a gene to be non-DE if it is not retained after filtering.).


```r
# Upregulated across most cell types.
up.de <- is.de > 0 & !is.na(is.de)
head(sort(rowMeans(up.de), decreasing=TRUE), 10)
```

```
## ENSG00000081059 ENSG00000137731 ENSG00000064886 ENSG00000066294 ENSG00000096060 
##       0.8181818       0.8181818       0.7272727       0.7272727       0.7272727 
## ENSG00000136161 ENSG00000158488 ENSG00000159674 ENSG00000169994 ENSG00000182866 
##       0.7272727       0.7272727       0.7272727       0.7272727       0.7272727
```


```r
# Downregulated across cell types.
down.de <- is.de < 0 & !is.na(is.de)
head(sort(rowMeans(down.de), decreasing=TRUE), 10)
```

```
## ENSG00000110492 ENSG00000019582 ENSG00000068079 ENSG00000100721 ENSG00000105369 
##       0.8181818       0.7272727       0.7272727       0.7272727       0.7272727 
## ENSG00000120833 ENSG00000162654 ENSG00000184489 ENSG00000196126 ENSG00000197872 
##       0.7272727       0.7272727       0.7272727       0.7272727       0.7272727
```

We further identify label-specific DE genes that are significant in our label of interest yet not DE in any other label. As hypothesis tests are not typically geared towards identifying genes that are not DE, we use an ad hoc approach where we consider a gene to be consistent with the null hypothesis for a label if it fails to be detected even at a generous FDR threshold of 50%.


```r
remotely.de <- decideTestsPerLabel(de.results, threshold=0.5)
not.de <- remotely.de==0 | is.na(remotely.de)

# first cluster in is.de
cx <- colnames(is.de)[1]
other.labels <- setdiff(colnames(not.de), cx)
unique.degs <- is.de[,cx]!=0 & rowMeans(not.de[,other.labels])==1
unique.degs <- names(which(unique.degs))
head(unique.degs)
```

```
## character(0)
```

```r
# 2nd cluster in is.de
cx <- colnames(is.de)[2]
other.labels <- setdiff(colnames(not.de), cx)
unique.degs <- is.de[,cx]!=0 & rowMeans(not.de[,other.labels])==1
unique.degs <- names(which(unique.degs))
```


```r
# Choosing the top-ranked gene for inspection:
de.inspec <- list()
de.inspec[[cx]] <- de.results[[cx]] 
de.inspec[[cx]] <- de.inspec[[cx]][order(de.inspec[[cx]]$PValue),]
de.inspec[[cx]] <- de.inspec[[cx]][rownames(de.inspec[[cx]]) %in% unique.degs,]

sizeFactors(summed.filt) <- NULL
plotExpression(logNormCounts(summed.filt), 
    features=rownames(de.inspec[[cx]])[1],
    x="source_name", colour_by="source_name", 
    other_fields="label") + 
    facet_wrap(~label)
```

<img src="multiSplComp_files/figure-html/unnamed-chunk-36-1.png" width="672" />

We also list the labels that were skipped due to the absence of replicates or contrasts. If it is necessary to extract statistics in the absence of replicates, several strategies can be applied such as reducing the complexity of the model or using a predefined value for the NB dispersion. We refer readers to the edgeR user’s guide for more details.


```r
print(metadata(de.results)$failed)
```

```
## [1] "c12" "c14" "c16" "c18" "c5"  "c7"  "c9"
```

## Differential abundance between conditions

### Overview

n a DA analysis, we test for significant changes in per-label cell abundance across conditions. This will reveal which cell types are depleted or enriched upon treatment, which is arguably just as interesting as changes in expression within each cell type. The DA analysis has a long history in flow cytometry (Finak et al. 2014; Lun, Richard, and Marioni 2017) where it is routinely used to examine the effects of different conditions on the composition of complex cell populations. By performing it here, we effectively treat scRNA-seq as a “super-FACS” technology for defining relevant subpopulations using the entire transcriptome.

We prepare for the DA analysis by quantifying the number of cells assigned to each label (or cluster).


```r
abundances <- table(merged$clusters.mnn, merged$Sample.Name2) 
abundances <- unclass(abundances) 
head(abundances)
```

```
##      
##       ETV6-RUNX1_1 ETV6-RUNX1_2 ETV6-RUNX1_3 ETV6-RUNX1_4 HHD_1 HHD_2 PBMMC_1
##   c1           937          722          165         1363   276   345      26
##   c10          103         2625          485          512   426  1042       9
##   c11           72          196          713           52    96   191      28
##   c12            0            0           11            7     3     1     396
##   c13            1            7           24           43    52    12     288
##   c14            3            3          135          253     1    15      64
##      
##       PBMMC_2 PBMMC_3 PRE-T_1 PRE-T_2
##   c1       27      41      37       2
##   c10      41     131     179     125
##   c11      18      64     171    1142
##   c12     152      13       2       2
##   c13     435     748       6     246
##   c14     611     142       2      96
```

Performing the DA analysis

Our DA analysis will again be performed with the edgeR package. This allows us to take advantage of the NB GLM methods to model overdispersed count data in the presence of limited replication - except that the counts are not of reads per gene, but of cells per label (Lun, Richard, and Marioni 2017). The aim is to share information across labels to improve our estimates of the biological variability in cell abundance between replicates.


```r
# Attaching some column metadata.
extra.info <- colData(merged)[match(colnames(abundances), merged$Sample.Name2),]
y.ab <- DGEList(abundances, samples=extra.info)
y.ab
```

```
## An object of class "DGEList"
## $counts
##      
##       ETV6-RUNX1_1 ETV6-RUNX1_2 ETV6-RUNX1_3 ETV6-RUNX1_4 HHD_1 HHD_2 PBMMC_1
##   c1           937          722          165         1363   276   345      26
##   c10          103         2625          485          512   426  1042       9
##   c11           72          196          713           52    96   191      28
##   c12            0            0           11            7     3     1     396
##   c13            1            7           24           43    52    12     288
##      
##       PBMMC_2 PBMMC_3 PRE-T_1 PRE-T_2
##   c1       27      41      37       2
##   c10      41     131     179     125
##   c11      18      64     171    1142
##   c12     152      13       2       2
##   c13     435     748       6     246
## 13 more rows ...
## 
## $samples
##              group lib.size norm.factors        batch        Run Sample.Name
## ETV6-RUNX1_1     1     2853            1 ETV6-RUNX1_1 SRR9264343  GSM3872434
## ETV6-RUNX1_2     1     6615            1 ETV6-RUNX1_2 SRR9264344  GSM3872435
## ETV6-RUNX1_3     1     4727            1 ETV6-RUNX1_3 SRR9264345  GSM3872436
## ETV6-RUNX1_4     1     5293            1 ETV6-RUNX1_4 SRR9264346  GSM3872437
## HHD_1            1     4551            1        HHD_1 SRR9264347  GSM3872438
##              source_name      block setName Sample.Name2 clusters.mnn
## ETV6-RUNX1_1  ETV6-RUNX1 ETV6-RUNX1   Caron ETV6-RUNX1_1           c3
## ETV6-RUNX1_2  ETV6-RUNX1 ETV6-RUNX1   Caron ETV6-RUNX1_2           c4
## ETV6-RUNX1_3  ETV6-RUNX1 ETV6-RUNX1   Caron ETV6-RUNX1_3          c17
## ETV6-RUNX1_4  ETV6-RUNX1 ETV6-RUNX1   Caron ETV6-RUNX1_4           c1
## HHD_1                HHD        HHD   Caron        HHD_1           c5
## 6 more rows ...
```

We filter out low-abundance labels as previously described. This avoids cluttering the result table with very rare subpopulations that contain only a handful of cells. For a DA analysis of cluster abundances, filtering is generally not required as most clusters will not be of low-abundance (otherwise there would not have been enough evidence to define the cluster in the first place).


```r
keep <- filterByExpr(y.ab, group=y.ab$samples$source_name)
y.ab <- y.ab[keep,]
summary(keep)
```

```
##    Mode    TRUE 
## logical      18
```

Unlike DE analyses, we do not perform an additional normalization step with calcNormFactors(). This means that we are only normalizing based on the “library size”, i.e., the total number of cells in each sample. Any changes we detect between conditions will subsequently represent differences in the proportion of cells in each cluster. The motivation behind this decision is discussed in more detail in Section 14.4.3.

Here, the log-fold change in our model refers to the change in cell abundance between sample groups, rather than the change in gene expression.


```r
design <- model.matrix(~factor(source_name), y.ab$samples)
```

We use the estimateDisp() function to estimate the NB dipersion for each cluster. We turn off the trend as we do not have enough points for its stable estimation.


```r
y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.9793  0.9793  0.9793  0.9793  0.9793  0.9793
```


```r
plotBCV(y.ab, cex=1)
```

<img src="multiSplComp_files/figure-html/unnamed-chunk-43-1.png" width="672" />

We repeat this process with the QL dispersion, again disabling the trend.


```r
fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    1.18    1.18    1.18    1.18    1.18    1.18
```


```r
summary(fit.ab$df.prior)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   120.9   212.5   212.5   199.0   212.5   212.5
```


```r
plotQLDisp(fit.ab, cex=1)
```

<img src="multiSplComp_files/figure-html/unnamed-chunk-46-1.png" width="672" />

We test for differences in abundance between sample groups using glmQLFTest().


```r
res <- glmQLFTest(fit.ab, coef=ncol(design))
summary(decideTests(res))
```

```
##        factor(source_name)PRE-T
## Down                          2
## NotSig                       16
## Up                            0
```


```r
topTags(res)
```

```
## Coefficient:  factor(source_name)PRE-T 
##         logFC   logCPM          F      PValue        FDR
## c5  -5.970062 16.35445 10.1640813 0.001806486 0.03251675
## c1  -5.153475 16.31480  8.3152532 0.004624676 0.04162208
## c13  3.227077 15.51131  5.9376070 0.016218199 0.09730919
## c7  -3.936424 15.34103  5.0729830 0.026032057 0.11714426
## c18  3.386502 11.42393  3.4451433 0.065774025 0.23608559
## c2   2.301146 16.85892  3.1424598 0.078695196 0.23608559
## c8   1.788105 15.51327  1.9509604 0.164938752 0.37206264
## c11  1.736525 15.89057  1.8117349 0.180717038 0.37206264
## c10 -1.963157 16.53987  1.7680054 0.186031318 0.37206264
## c16  1.144352 11.89816  0.7203937 0.397624426 0.71572397
```

### Handling composition effects

#### Background

As mentioned above, we do not use calcNormFactors() in our default DA analysis. This normalization step assumes that most of the input features are not different between conditions. While this assumption is reasonable for most types of gene expression data, it is generally too strong for cell type abundance - most experiments consist of only a few cell types that may all change in abundance upon perturbation. Thus, our default approach is to only normalize based on the total number of cells in each sample, which means that we are effectively testing for differential proportions between conditions.

Unfortunately, the use of the total number of cells leaves us susceptible to composition effects. For example, a large increase in abundance for one cell subpopulation will introduce decreases in proportion for all other subpopulations - which is technically correct, but may be misleading if one concludes that those other subpopulations are decreasing in abundance of their own volition. If composition biases are proving problematic for interpretation of DA results, we have several avenues for removing them or mitigating their impact by leveraging a priori biological knowledge.
14.4.3.2 Assuming most labels do not change

If it is possible to assume that most labels (i.e., cell types) do not change in abundance, we can use calcNormFactors() to compute normalization factors.


```r
y.ab2 <- calcNormFactors(y.ab)
y.ab2$samples$norm.factors
```

```
##  [1] 0.6775040 0.7238884 0.9769490 1.1290651 0.9268008 0.9656680 1.2257770
##  [8] 1.4773071 1.4480029 0.7174753 1.0978745
```

We then proceed with the remainder of the edgeR analysis, shown below in condensed format. A shift of positive log-fold changes towards zero is consistent with the removal of composition biases.


```r
y.ab2 <- estimateDisp(y.ab2, design, trend="none")
fit.ab2 <- glmQLFit(y.ab2, design, robust=TRUE, abundance.trend=FALSE)
res2 <- glmQLFTest(fit.ab2, coef=ncol(design))
topTags(res2, n=10)
```

```
## Coefficient:  factor(source_name)PRE-T 
##         logFC   logCPM          F      PValue        FDR
## c5  -6.274363 16.53184 10.3342094 0.001658929 0.02986073
## c1  -5.017465 16.56753  7.6891289 0.006399471 0.05759524
## c13  3.118403 15.14410  5.4457475 0.021199590 0.11692701
## c7  -3.826268 15.01359  5.0763350 0.025983780 0.11692701
## c18  3.430538 11.09678  3.0342207 0.083965998 0.30227759
## c2   2.072243 17.09976  2.4964424 0.116610952 0.34983286
## c8   1.811393 15.79196  1.9041067 0.170063044 0.40510432
## c10 -2.040034 16.76314  1.8173604 0.180046364 0.40510432
## c11  1.530994 15.91093  1.3570098 0.246257110 0.49251422
## c17 -0.960683 15.18727  0.4518357 0.502695247 0.83541716
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
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] edgeR_3.32.1                limma_3.46.0               
##  [3] batchelor_1.6.3             Cairo_1.5-12.2             
##  [5] DT_0.18                     dplyr_1.0.5                
##  [7] scran_1.18.7                scater_1.18.6              
##  [9] SingleCellExperiment_1.12.0 SummarizedExperiment_1.20.0
## [11] Biobase_2.50.0              GenomicRanges_1.42.0       
## [13] GenomeInfoDb_1.26.7         IRanges_2.24.1             
## [15] S4Vectors_0.28.1            BiocGenerics_0.36.1        
## [17] MatrixGenerics_1.2.1        matrixStats_0.58.0         
## [19] ggplot2_3.3.3               knitr_1.32                 
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-7              RColorBrewer_1.1-2       
##  [3] tools_4.0.3               bslib_0.2.4              
##  [5] utf8_1.2.1                R6_2.5.0                 
##  [7] irlba_2.3.3               ResidualMatrix_1.0.0     
##  [9] vipor_0.4.5               uwot_0.1.10              
## [11] DBI_1.1.1                 colorspace_2.0-0         
## [13] withr_2.4.2               tidyselect_1.1.1         
## [15] gridExtra_2.3             compiler_4.0.3           
## [17] cli_2.4.0                 BiocNeighbors_1.8.2      
## [19] DelayedArray_0.16.3       labeling_0.4.2           
## [21] bookdown_0.22             sass_0.3.1               
## [23] scales_1.1.1              stringr_1.4.0            
## [25] digest_0.6.27             rmarkdown_2.7            
## [27] XVector_0.30.0            pkgconfig_2.0.3          
## [29] htmltools_0.5.1.1         sparseMatrixStats_1.2.1  
## [31] fastmap_1.1.0             highr_0.9                
## [33] htmlwidgets_1.5.3         rlang_0.4.10             
## [35] rstudioapi_0.13           shiny_1.6.0              
## [37] DelayedMatrixStats_1.12.3 farver_2.1.0             
## [39] jquerylib_0.1.3           generics_0.1.0           
## [41] jsonlite_1.7.2            crosstalk_1.1.1          
## [43] BiocParallel_1.24.1       RCurl_1.98-1.3           
## [45] magrittr_2.0.1            BiocSingular_1.6.0       
## [47] GenomeInfoDbData_1.2.4    scuttle_1.0.4            
## [49] Matrix_1.3-2              Rcpp_1.0.6               
## [51] ggbeeswarm_0.6.0          munsell_0.5.0            
## [53] fansi_0.4.2               viridis_0.6.0            
## [55] lifecycle_1.0.0           stringi_1.5.3            
## [57] yaml_2.2.1                zlibbioc_1.36.0          
## [59] Rtsne_0.15                grid_4.0.3               
## [61] promises_1.2.0.1          dqrng_0.3.0              
## [63] crayon_1.4.1              lattice_0.20-44          
## [65] splines_4.0.3             cowplot_1.1.1            
## [67] beachmat_2.6.4            locfit_1.5-9.4           
## [69] pillar_1.6.0              igraph_1.2.6             
## [71] codetools_0.2-18          glue_1.4.2               
## [73] evaluate_0.14             httpuv_1.5.5             
## [75] vctrs_0.3.7               gtable_0.3.0             
## [77] purrr_0.3.4               assertthat_0.2.1         
## [79] xfun_0.22                 mime_0.10                
## [81] rsvd_1.0.5                xtable_1.8-4             
## [83] RSpectra_0.16-0           later_1.2.0              
## [85] viridisLite_0.4.0         pheatmap_1.0.12          
## [87] tibble_3.1.1              beeswarm_0.3.1           
## [89] bluster_1.0.0             statmod_1.4.35           
## [91] ellipsis_0.3.2
```
