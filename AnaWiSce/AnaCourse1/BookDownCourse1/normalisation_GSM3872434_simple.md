---
title: "CRUK CI Summer School 2020 - introduction to single-cell RNA-seq analysis"
subtitle: 'Normalisation - GSM3872434 sample'
author: "Stephane Ballereau, Zeynep Kalender Atak, Katarzyna Kania"
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
  splSetToGet: "dummy"
  setName: "GSM3872434"
  setSuf: "_allCells"
  dsiSuf: 'dummy'
---

# Normalisation - GSM3872434 set {#NormalisationGSM3872434Top}




Sources: chapters on normalisation in the [OSCA book](https://osca.bioconductor.org/normalization.html) and the ['Hemberg group material'](https://scrnaseq-course.cog.sanger.ac.uk/website/index.html).

Why normalise?

Systematic differences in sequencing coverage between libraries occur because of low input material, differences in cDNA capture and PCR amplification.
Normalisation removes such differences so that differences between cells are not technical but biological, allowing meaningful comparison of expression profiles between cells. Normalisation and batch correction have different aims. Normalisation addresses technical differences only, while batch correction considers both technical and biological differences.


```r
normPlotDirBit <- "Plots/Norm"

projDir <- params$projDir
dirRel <- params$dirRel
outDirBit <- params$outDirBit
setName <- params$setName
setSuf <- params$setSuf

if(params$bookType == "mk") {
	setName <- "GSM3872434"
	setSuf <- "_allCells"
	dirRel <- ".."
}
dir.create(sprintf("%s/%s/%s", projDir, outDirBit, normPlotDirBit),
	   showWarnings = FALSE)
```



Load object 


```r
# Read object in:
tmpFn <- sprintf("%s/%s/Robjects/%s_sce_nz_postQc%s.Rds",
		 projDir, outDirBit, "caron", setSuf)
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
## assays(1): counts
## rownames(18431): ENSG00000238009 ENSG00000237491 ... ENSG00000275063
##   ENSG00000271254
## rowData names(11): ensembl_gene_id external_gene_name ... detected
##   gene_sparsity
## colnames: NULL
## colData names(16): Sample Barcode ... outlier cell_sparsity
## reducedDimNames(0):
## altExpNames(0):
```

Select cells for GSM3872434:


```r
##setSuf <- "_5hCellPerSpl"
##nbCells <- 500
#setSuf <- "_1kCellPerSpl"
#nbCells <- 1000
setSuf <- "_GSM3872434"
##nbCells <- 500
```


```r
# have new list of cell barcodes for each sample
sce.nz.master <- sce
vec.bc <- colData(sce.nz.master) %>%
	data.frame() %>%
	filter(Sample.Name == "GSM3872434") %>%
	group_by(Sample.Name) %>%
	##sample_n(nbCells) %>%
	pull(Barcode)
```

Number of cells in the sample:


```r
table(colData(sce.nz.master)$Barcode %in% vec.bc)
```

```
## 
## FALSE  TRUE 
## 44977  2853
```

Subset cells from the SCE object:


```r
tmpInd <- which(colData(sce.nz.master)$Barcode %in% vec.bc)
sce <- sce.nz.master[,tmpInd]
sce
```

```
## class: SingleCellExperiment 
## dim: 18431 2853 
## metadata(0):
## assays(1): counts
## rownames(18431): ENSG00000238009 ENSG00000237491 ... ENSG00000275063
##   ENSG00000271254
## rowData names(11): ensembl_gene_id external_gene_name ... detected
##   gene_sparsity
## colnames: NULL
## colData names(16): Sample Barcode ... outlier cell_sparsity
## reducedDimNames(0):
## altExpNames(0):
```

Check columns data:


```r
head(colData(sce))
```

```
## DataFrame with 6 rows and 16 columns
##                   Sample            Barcode         Run Sample.Name source_name
##              <character>        <character> <character> <character>    <factor>
## 1 /ssd/personal/baller.. AAACCTGAGACTTTCG-1  SRR9264343  GSM3872434  ETV6-RUNX1
## 2 /ssd/personal/baller.. AAACCTGGTCTTCAAG-1  SRR9264343  GSM3872434  ETV6-RUNX1
## 3 /ssd/personal/baller.. AAACCTGGTGTTGAGG-1  SRR9264343  GSM3872434  ETV6-RUNX1
## 4 /ssd/personal/baller.. AAACCTGTCCCAAGTA-1  SRR9264343  GSM3872434  ETV6-RUNX1
## 5 /ssd/personal/baller.. AAACCTGTCGAATGCT-1  SRR9264343  GSM3872434  ETV6-RUNX1
## 6 /ssd/personal/baller.. AAACGGGCACCATCCT-1  SRR9264343  GSM3872434  ETV6-RUNX1
##         sum  detected subsets_Mito_sum subsets_Mito_detected
##   <numeric> <integer>        <numeric>             <integer>
## 1      6462      1996              290                    12
## 2     11706      3079              567                    12
## 3      7981      2511              423                    12
## 4      8354      2307              519                    13
## 5      1373       701               91                    11
## 6      2142       985              102                    10
##   subsets_Mito_percent     total      block     setName   discard   outlier
##              <numeric> <numeric>   <factor> <character> <logical> <logical>
## 1              4.48777      6462 ETV6-RUNX1       Caron     FALSE     FALSE
## 2              4.84367     11706 ETV6-RUNX1       Caron     FALSE     FALSE
## 3              5.30009      7981 ETV6-RUNX1       Caron     FALSE     FALSE
## 4              6.21259      8354 ETV6-RUNX1       Caron     FALSE     FALSE
## 5              6.62782      1373 ETV6-RUNX1       Caron     FALSE     FALSE
## 6              4.76190      2142 ETV6-RUNX1       Caron     FALSE     FALSE
##   cell_sparsity
##       <numeric>
## 1      0.907503
## 2      0.857361
## 3      0.883637
## 4      0.893091
## 5      0.967515
## 6      0.954400
```

```r
table(colData(sce)$Sample.Name)
```

```
## 
## GSM3872434 
##       2853
```

<!--
# mind that genes were filtered using all cells, not just those sampled here.
-->

We write the R object to GSM3872434_sce_nz_postQc_GSM3872434.Rds.


```r
# Write object to file
tmpFn <- sprintf("%s/%s/Robjects/%s_sce_nz_postQc%s.Rds",
		 projDir, outDirBit, setName, setSuf)
saveRDS(sce, tmpFn)
```


```r
# Write object to file
tmpFn <- sprintf("%s/%s/Robjects/%s_sce_nz_postQc%s.Rds",
		 projDir, outDirBit, setName, setSuf)
sce <- readRDS(tmpFn)
```

## Scaling normalization

In scaling normalization, the “normalization factor” is an estimate of the library size relative to the other cells. Steps usually include: computation of a cell-specific 'scaling' or 'size' factor that represents the relative bias in that cell and division of all counts for the cell by that factor to remove that bias. Assumption: any cell specific bias will affect genes the same way.

Scaling methods typically generate normalised counts-per-million (CPM) or transcripts-per-million (TPM) values that address the effect of sequencing depth. These values however typically have a variance that changes with their mean (heteroscedasticity) while most statistical methods assume a stable variance, which does not vary with the mean (homoscedasticity). A widely used 'variance stabilising transformation' is the log transformation (often log2). This works fine for highly expressed genes (as in bulk RNA-seq) but less so for sparse scRNA-seq data.

### CPM

Convert raw counts to counts-per-million (CPM) for each cell by dividing counts by the library size then multiplying by 1.000.000. Mind that this does not adress compositional bias caused by highly expressed genes that are also differentially expressed between cells. In `scater` CPMs are computed with the following code:


```r
calc_cpm <- function (expr_mat, spikes = NULL) 
{
    norm_factor <- colSums(expr_mat[-spikes, ])
    return(t(t(expr_mat)/norm_factor)) * 10^6
}
```

We will use `scater`'s `calculateCPM()`.

### DESeq's size factor

For each gene, compute geometric mean across cells. For each cell compute for each gene the ratio of its expression to its geometric mean, and derive the cell's size factor as the median ratio across genes. Not suitable for sparse scRNA-seq data as the geometric is computed on non-zero values only. This method is also known as 'Relative Log Expression' (RLE) in `edgeR` and `scater`. 

Example code:


```r
calc_sf <- function (expr_mat, spikes = NULL) 
{
    geomeans <- exp(rowMeans(log(expr_mat[-spikes, ])))
    SF <- function(cnts) {
        median((cnts/geomeans)[(is.finite(geomeans) &
				geomeans > 0)])
    }
    norm_factor <- apply(expr_mat[-spikes, ], 2, SF)
    return(t(t(expr_mat)/norm_factor))
}
```

### Weighted Trimmed mean of M-values

To compute weighted Trimmed mean of M-values (TMM), a given cell is chosen as a reference to use in computation for other cells. The M-values are gene-wise log2-fold changes between cells. Trimming entails the removal of the top and bottom 30% of values. The size factor is computed as the average for the remaining cells with a weight according to inverse variances. This method assumes that most genes are not differentially expressed, and the 40% of genes left after trimming may include many zero counts.


```r
sizeFactors(sce) <- edgeR::calcNormFactors(counts(sce), method = "TMM")
```

### Library size factor distribution

For each cell, the library size factor is proportional to the library size such that the average size factor across cells is one.

Advantage: normalised counts are on the same scale as the initial counts.

Compute size factors:


```r
lib.sf <- librarySizeFactors(sce)
summary(lib.sf)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.09666 0.66197 0.92591 1.00000 1.23957 5.49973
```

Size factor distribution: wide range, typical of scRNA-seq data.


```r
hist(log10(lib.sf), xlab="Log10[Size factor]", col='grey80')
```

<img src="normalisation_GSM3872434_simple_files/figure-html/librarySizeFactors_hist_norm_GSM3872434_allCells-1.png" width="672" />

Assumption: absence of compositional bias; differential expression between two cells is balanced: upregulation in some genes is accompanied by downregulation of other genes. Not observed.

Inaccurate normalisation due to unaccounted-for composition bias affects the size of the log fold change measured between clusters, but less so the clustering itself. It is thus sufficient to identify clusters and top marker genes.

### Deconvolution

Composition bias occurs when differential expression beteween two samples or here cells is not balanced. For a fixed library size, identical in both cells, upregulation of one gene in a cell will means fewer UMIs can be assigned to other genes, which would then appear down regulated. Even if library sizes are allowed to differ, with that for the cell with upregulation being higher, scaling normalisation will reduce normalised counts. Non-upregulated would therefore also appear downregulated. 

For bulk RNA-seq, composition bias is removed by assuming that most genes are not differentially expressed between samples, so that differences in non-DE genes would amount to the bias, and used to compute size factors.

Given the sparsity of scRNA-seq data, the methods are not appropriate.

The method below increases read counts by pooling cells into groups, computing size factors within each of these groups and scaling them so they are comparable across clusters. This process is repeated many times, changing pools each time to collect several size factors for each cell, from which is derived a single value for that cell.

<!--
see DESeq2 estimateSizeFactorsFromMatrix
see edgeR calcNormFactors
-->


```r
tmpFn <- sprintf("%s/Images/scran_Fig3.png", "..")
knitr::include_graphics(tmpFn, auto_pdf = TRUE)
```

<img src="../Images/scran_Fig3.png" width="398" />

```r
rm(tmpFn)
```

Cluster cells then normalise.

#### Cluster cells


```r
set.seed(100) # clusters with PCA from irlba with approximation
clust <- quickCluster(sce) # slow with all cells.

# write to file
tmpFn <- sprintf("%s/%s/Robjects/%s_sce_nz_quickClus%s.Rds",
		 projDir, outDirBit, setName, setSuf)
saveRDS(clust, tmpFn)
```


```r
# read from file
tmpFn <- sprintf("%s/%s/Robjects/%s_sce_nz_quickClus%s.Rds",
		 projDir, outDirBit, setName, setSuf)
clust <- readRDS(tmpFn)
table(clust)
```

```
## clust
##   1   2   3   4   5   6   7   8   9 
## 232 193 617 315 407 380 398 189 122
```

#### Compute size factors


```r
#deconv.sf <- calculateSumFactors(sce, cluster=clust)
sce <- computeSumFactors(sce, cluster=clust, min.mean=0.1)
deconv.sf <- sizeFactors(sce)

# write to file
tmpFn <- sprintf("%s/%s/Robjects/%s_sce_nz_deconvSf%s.Rds",
		 projDir, outDirBit, setName, setSuf)
saveRDS(deconv.sf, tmpFn)
```


```r
# read from file
tmpFn <- sprintf("%s/%s/Robjects/%s_sce_nz_deconvSf%s.Rds",
		 projDir, outDirBit, setName, setSuf)
deconv.sf <- readRDS(tmpFn)

summary(deconv.sf)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.08207 0.67672 0.91982 1.00000 1.22628 4.06014
```

Plot deconvolution size factors against library size factors:




```r
deconvDf <- data.frame(lib.sf, deconv.sf,
			"source_name" = sce$source_name,
			"sum" = sce$sum,
			"mito_content" = sce$subsets_Mito_percent,
			"cell_sparsity" = sce$cell_sparsity)
```


```r
sp <- ggplot(deconvDf, aes(x=lib.sf, y=deconv.sf, col=source_name)) +
  geom_point()

# Split by sample type:
#sp + facet_wrap(~source_name)
```




```r
sp <- ggplot(deconvDf, aes(x=lib.sf, y=deconv.sf, col=cell_sparsity)) +
  geom_point()
sp
```

<img src="normalisation_GSM3872434_simple_files/figure-html/scatter_deconvSf_libSf_colBy_cellSpars_norm_GSM3872434_allCells-1.png" width="672" />

#### Apply size factors

For each cell, raw counts for genes are divided by the size factor for that cell and log-transformed so downstream analyses focus on genes with strong relative differences. We use `scater::logNormCounts()`.


```r
sce <- logNormCounts(sce) # adds logcounts
print(assayNames(sce))
```

```
## [1] "counts"    "logcounts"
```

#### Save object


```r
# write to file
tmpFn <- sprintf("%s/%s/Robjects/%s_sce_nz_postDeconv%s.Rds", projDir, outDirBit, setName, setSuf)
saveRDS(sce, tmpFn)
```


```r
sceDeconv <- sce
```

## SCTransform

<!--
https://rawgit.com/ChristophH/sctransform/master/inst/doc/variance_stabilizing_transformation.html

vars.to.regress = c("S.Score", "G2M.Score")
vars.to.regress = c("percentMito","Sex")
-->

With scaling normalisation a correlation remains between the mean and variation of expression (heteroskedasticity). This affects downstream dimensionality reduction as the few main new dimensions are usually correlated with library size. SCTransform addresses the issue by regressing library size out of raw counts and providing residuals to use as normalized and variance-stabilized expression values in downstream analysis. We will use the [sctransform vignette](https://cran.r-project.org/web/packages/sctransform/index.html).


```r
counts <- counts(sce)
print(class(counts))
```

```
## [1] "dgCMatrix"
## attr(,"package")
## [1] "Matrix"
```

```r
print(dim(counts))
```

```
## [1] 18431  2853
```

```r
colnames(counts) <- colData(sce)$Barcode
```

### Inspect data

We will now calculate some properties and visually inspect the data. Our main interest is in the general trends not in individual outliers. Neither genes nor cells that stand out are important at this step, but we focus on the global trends.

Derive gene and cell attributes from the UMI matrix.


```r
gene_attr <- data.frame(mean = rowMeans(counts), 
                        detection_rate = rowMeans(counts > 0),
                        var = apply(counts, 1, var))
gene_attr$log_mean <- log10(gene_attr$mean)
gene_attr$log_var <- log10(gene_attr$var)
rownames(gene_attr) <- rownames(counts)
cell_attr <- data.frame(n_umi = colSums(counts),
                        n_gene = colSums(counts > 0))
rownames(cell_attr) <- colnames(counts)
```


```r
dim(gene_attr)
```

```
## [1] 18431     5
```

```r
head(gene_attr)
```

```
##                         mean detection_rate          var   log_mean    log_var
## ENSG00000238009 0.0007010165   0.0007010165 0.0007007707 -3.1542718 -3.1544241
## ENSG00000237491 0.0364528566   0.0357518402 0.0365388860 -1.4382684 -1.4372447
## ENSG00000225880 0.0318962496   0.0304942166 0.0336947550 -1.4962604 -1.4724377
## ENSG00000230368 0.0140203295   0.0136698212 0.0145298692 -1.8532418 -1.8377383
## ENSG00000230699 0.0021030494   0.0021030494 0.0020993624 -2.6771505 -2.6779126
## ENSG00000188976 0.2958289520   0.2509638977 0.3065639427 -0.5289593 -0.5134789
```


```r
dim(cell_attr)
```

```
## [1] 2853    2
```

```r
head(cell_attr)
```

```
##                    n_umi n_gene
## AAACCTGAGACTTTCG-1  6462   1996
## AAACCTGGTCTTCAAG-1 11705   3078
## AAACCTGGTGTTGAGG-1  7981   2511
## AAACCTGTCCCAAGTA-1  8354   2307
## AAACCTGTCGAATGCT-1  1373    701
## AAACGGGCACCATCCT-1  2141    984
```

Mean-variance relationship

For the genes, we can see that up to a mean UMI count of 0 the variance follows the line through the origin with slop one, i.e. variance and mean are roughly equal as expected under a Poisson model. However, genes with a higher average UMI count show overdispersion compared to Poisson.


```r
ggplot(gene_attr, aes(log_mean, log_var)) + 
  geom_point(alpha=0.3, shape=16) + 
  geom_density_2d(size = 0.3) +
  geom_abline(intercept = 0, slope = 1, color='red')
```

<img src="normalisation_GSM3872434_simple_files/figure-html/attr_plot_sct_GSM3872434_allCells-1.png" width="672" />

Mean-detection-rate relationship

In line with the previous plot, we see a lower than expected detection rate in the medium expression range. However, for the highly expressed genes, the rate is at or very close to 1.0 suggesting that there is no zero-inflation in the counts for those genes and that zero-inflation is a result of overdispersion, rather than an independent systematic bias.


```r
# add the expected detection rate under Poisson model
x = seq(from = -3, to = 2, length.out = 1000)
poisson_model <- data.frame(log_mean = x, detection_rate = 1 - dpois(0, lambda = 10^x))
ggplot(gene_attr, aes(log_mean, detection_rate)) + 
  geom_point(alpha=0.3, shape=16) + 
  geom_line(data=poisson_model, color='red') +
  theme_gray(base_size = 8)
```

<img src="normalisation_GSM3872434_simple_files/figure-html/scatter_detecRate_logMean_sct_GSM3872434_allCells-1.png" width="672" />




### Transformation

"Based on the observations above, which are not unique to this particular data set, we propose to model the expression of each gene as a negative binomial random variable with a mean that depends on other variables. Here the other variables can be used to model the differences in sequencing depth between cells and are used as independent variables in a regression model. In order to avoid overfitting, we will first fit model parameters per gene, and then use the relationship between gene mean and parameter values to fit parameters, thereby combining information across genes. Given the fitted model parameters, we transform each observed UMI count into a Pearson residual which can be interpreted as the number of standard deviations an observed count was away from the expected mean. If the model accurately describes the mean-variance relationship and the dependency of mean and latent factors, then the result should have mean zero and a stable variance across the range of expression." [sctransform vignette](https://cran.r-project.org/web/packages/sctransform/index.html).

Estimate model parameters and transform data

The vst function estimates model parameters and performs the variance stabilizing transformation. Here we use the log10 of the total UMI counts of a cell as variable for sequencing depth for each cell. After data transformation we plot the model parameters as a function of gene mean (geometric mean).


```r
print(dim(counts))
```

```
## [1] 18431  2853
```

```r
# We use the Future API for parallel processing; set parameters here
future::plan(strategy = 'multicore', workers = 4)
options(future.globals.maxSize = 10 * 1024 ^ 3)

set.seed(44)
vst_out <- sctransform::vst(counts,
                            latent_var = c('log_umi'),
                            return_gene_attr = TRUE,
                            return_cell_attr = TRUE,
                            show_progress = FALSE)
sctransform::plot_model_pars(vst_out)
```

<img src="normalisation_GSM3872434_simple_files/figure-html/comp_sct_GSM3872434_allCells-1.png" width="672" />

Inspect model:


```r
print(vst_out$model_str)
```

```
## [1] "y ~ log_umi"
```

We will look at several genes in more detail.


```r
#sctransform::plot_model(vst_out, counts, c('MALAT1', 'RPL10', 'FTL'), plot_residual = TRUE)
rowData(sce) %>%
  as.data.frame %>%
  filter(Symbol %in% c('MALAT1', 'RPL10', 'FTL'))
```

```
##                 ensembl_gene_id external_gene_name chromosome_name
## ENSG00000147403 ENSG00000147403              RPL10               X
## ENSG00000251562 ENSG00000251562             MALAT1              11
## ENSG00000087086 ENSG00000087086                FTL              19
##                 start_position end_position strand Symbol            Type
## ENSG00000147403      154389955    154409168      1  RPL10 Gene Expression
## ENSG00000251562       65497688     65506516      1 MALAT1 Gene Expression
## ENSG00000087086       48965309     48966879      1    FTL Gene Expression
##                      mean detected gene_sparsity
## ENSG00000147403  48.85548 99.28030   0.011540874
## ENSG00000251562 193.05529 99.09786   0.005352289
## ENSG00000087086  15.48621 96.42481   0.085532093
```

```r
sctransform::plot_model(vst_out,
                        counts,
                        c('ENSG00000251562', 'ENSG00000147403', 'ENSG00000087086'),
                        plot_residual = TRUE)
```

<img src="normalisation_GSM3872434_simple_files/figure-html/plot_model_1_sct_GSM3872434_allCells-1.png" width="672" />


```r
sctransform::plot_model(vst_out,
                        counts,
                        c('ENSG00000087086'),
                        plot_residual = TRUE,
                        show_nr = TRUE,
                        arrange_vertical = FALSE)
```

<img src="normalisation_GSM3872434_simple_files/figure-html/plot_model_2_sct_GSM3872434_allCells-1.png" width="672" />

Distribution of residual mean:


```r
ggplot(vst_out$gene_attr, aes(residual_mean)) +
  geom_histogram(binwidth=0.01)
```

<img src="normalisation_GSM3872434_simple_files/figure-html/plot_model_resMean_sct_GSM3872434_allCells-1.png" width="672" />

Distribution of residual variance:


```r
ggplot(vst_out$gene_attr, aes(residual_variance)) +
  geom_histogram(binwidth=0.1) +
  geom_vline(xintercept=1, color='red') +
  xlim(0, 10)
```

<img src="normalisation_GSM3872434_simple_files/figure-html/plot_model_resVar_sct_GSM3872434_allCells-1.png" width="672" />

Variance against mean (residuals):


```r
ggplot(vst_out$gene_attr, aes(x=residual_mean, y=residual_variance)) +
	geom_point(alpha=0.3, shape=16) + 
	xlim(0, 2.5) +
	ylim(0, 10) +
	geom_density_2d()
```

<img src="normalisation_GSM3872434_simple_files/figure-html/plot_model_resVar_resMean_sct_GSM3872434_allCells-1.png" width="672" />

Variance against mean (genes):


```r
ggplot(vst_out$gene_attr, aes(log10(gmean), residual_variance)) +
	geom_point(alpha=0.3, shape=16) +
  geom_density_2d(size = 0.3)
```

<img src="normalisation_GSM3872434_simple_files/figure-html/plot_model_resVar_gMean_sct_GSM3872434_allCells-1.png" width="672" />

Check genes with large residual variance:


```r
dd <- vst_out$gene_attr %>%
	arrange(-residual_variance) %>%
	slice_head(n = 22) %>%
	mutate(across(where(is.numeric), round, 2))

dd %>% tibble::rownames_to_column("ensembl_gene_id") %>%
	left_join(as.data.frame(rowData(sce))[,c("ensembl_gene_id", "Symbol")], "ensembl_gene_id")
```

```
##    ensembl_gene_id detection_rate gmean variance residual_mean
## 1  ENSG00000197061           0.72  2.04   534.32          1.90
## 2  ENSG00000175063           0.21  0.36    26.18          0.68
## 3  ENSG00000164104           0.66  1.74   120.31          1.06
## 4  ENSG00000100097           0.06  0.07     3.85          0.31
## 5  ENSG00000008517           0.01  0.02     0.87          0.20
## 6  ENSG00000123416           0.91  5.22   447.36          1.07
## 7  ENSG00000211679           0.08  0.09     1.46          0.30
## 8  ENSG00000102970           0.03  0.03     2.22          0.15
## 9  ENSG00000170540           0.54  0.86    38.67          0.47
## 10 ENSG00000164611           0.21  0.31    12.05          0.41
## 11 ENSG00000128322           0.33  0.61    11.86          0.69
## 12 ENSG00000244734           0.45  0.49 85356.97          0.15
## 13 ENSG00000188536           0.16  0.15  8088.98          0.11
## 14 ENSG00000206172           0.12  0.11  4995.50          0.09
## 15 ENSG00000271503           0.02  0.02     0.41          0.12
## 16 ENSG00000131747           0.16  0.27     7.92          0.34
## 17 ENSG00000187837           0.67  1.53    21.13          0.78
## 18 ENSG00000026025           0.69  1.92    58.23          0.73
## 19 ENSG00000115523           0.00  0.01     0.24          0.06
## 20 ENSG00000145649           0.00  0.01     0.06          0.08
## 21 ENSG00000158578           0.01  0.01     1.67          0.06
## 22 ENSG00000117399           0.11  0.16     4.07          0.25
##    residual_variance   Symbol
## 1              65.39 HIST1H4C
## 2              20.67    UBE2C
## 3              16.78    HMGB2
## 4              16.15   LGALS1
## 5              15.36     IL32
## 6              13.55   TUBA1B
## 7              12.64    IGLC3
## 8              10.97    CCL17
## 9              10.80  ARL6IP1
## 10             10.77    PTTG1
## 11             10.17    IGLL1
## 12              9.95      HBB
## 13              9.83     HBA2
## 14              9.75     HBA1
## 15              9.02     CCL5
## 16              8.63    TOP2A
## 17              8.44 HIST1H1C
## 18              8.42      VIM
## 19              7.85     GNLY
## 20              7.79     GZMA
## 21              7.78    ALAS2
## 22              7.53    CDC20
```

Write outcome to file:


```r
# write to file
tmpFn <- sprintf("%s/%s/Robjects/%s_sce_nz_vst_out%s.Rds",
		 projDir, outDirBit, setName, setSuf)
saveRDS(vst_out, tmpFn)
```

Check transformed values.


```r
print(dim(vst_out$y))
```

```
## [1] 14214  2853
```

```r
vst_out$y[1:10,1:5]
```

```
##                 AAACCTGAGACTTTCG-1 AAACCTGGTCTTCAAG-1 AAACCTGGTGTTGAGG-1
## ENSG00000237491        -0.20121617        -0.25944999        -0.22046915
## ENSG00000225880        -0.18668543        -0.24052717        -0.20448848
## ENSG00000230368        -0.12439501        -0.16164785        -0.13667117
## ENSG00000230699        -0.04842402        -0.06042004        -0.05240396
## ENSG00000188976         1.03374030         1.56177764         0.80159808
## ENSG00000187961        -0.05363629        -0.06698983        -0.05806420
## ENSG00000188290        -0.09174430        -0.12027862        -0.10110587
## ENSG00000187608         2.45936165        -0.77304781        -0.65435652
## ENSG00000188157        -0.07078407        -0.09039073        -0.07724886
## ENSG00000131591        -0.16874897        -0.21715523        -0.18475767
##                 AAACCTGTCCCAAGTA-1 AAACCTGTCGAATGCT-1
## ENSG00000237491        -0.22484005        -0.10084504
## ENSG00000225880        -0.20852982        -0.09378858
## ENSG00000230368        -0.13946384        -0.06115474
## ENSG00000230699        -0.05330481        -0.02694137
## ENSG00000188976         0.75303608        -0.27198451
## ENSG00000187961        -0.05906678        -0.02977692
## ENSG00000188290        -0.10324123        -0.04420393
## ENSG00000187608        -0.66778017        -0.28420778
## ENSG00000188157        -0.07871818        -0.03690960
## ENSG00000131591        -0.18839118        -0.08509713
```

Genes that are expressed in fewer than 5 cells are not used and not returned,
so to add vst_out$y as an assay we need to remove the missing genes.


```r
# https://github.com/ChristophH/sctransform/issues/27
sceOrig <- sce
sceOrig
```

```
## class: SingleCellExperiment 
## dim: 18431 2853 
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
tmpInd <- which(rownames(sce) %in% rownames(vst_out$y))
cols.meta <- colData(sceOrig)
rows.meta <- rowData(sceOrig)
new.counts <- counts(sceOrig)[tmpInd, ]
sce <- SingleCellExperiment(list(counts=new.counts))

# reset the column data on the new object
colData(sce) <- cols.meta
rowData(sce) <- rows.meta[tmpInd, ]
assayNames(sce)
```

```
## [1] "counts"
```


```r
sce
```

```
## class: SingleCellExperiment 
## dim: 14214 2853 
## metadata(0):
## assays(1): counts
## rownames(14214): ENSG00000237491 ENSG00000225880 ... ENSG00000278817
##   ENSG00000271254
## rowData names(11): ensembl_gene_id external_gene_name ... detected
##   gene_sparsity
## colnames: NULL
## colData names(17): Sample Barcode ... cell_sparsity sizeFactor
## reducedDimNames(0):
## altExpNames(0):
```

```r
vstMat <- as(vst_out$y[rownames(sce),], "dgCMatrix")
all(colnames(vstMat) == sce$Barcode)
```

```
## [1] TRUE
```

```r
colnames(vstMat) <- NULL
assay(sce, "sctrans_norm") <- vstMat # as(vst_out$y[rownames(sce),], "dgCMatrix")
#assayNames(sce)
```

### Save SCE object


```r
# write to file
tmpFn <- sprintf("%s/%s/Robjects/%s_sce_nz_postSct%s.Rds",
		 projDir, outDirBit, setName, setSuf)
saveRDS(sce, tmpFn)
```

## Visualisation

### log raw counts


```r
typeNorm <- "logRaw"

#setSuf <- "_5kCellPerSpl"
options(BiocSingularParam.default=IrlbaParam())

assay(sce, "logcounts_raw") <- log2(counts(sce) + 1)

tmp <- runPCA(
  sce[,],
  exprs_values = "logcounts_raw"
)
```

PCA plot for the logRaw counts in the GSM3872434 set.


```r
p <- plotPCA(
    tmp,
    colour_by = "Sample.Name",
    size_by = "sum",
    shape_by = "source_name"
) + ggtitle(sprintf("PCA plot: %s", typeNorm))

# write plot to file:
tmpFn <- sprintf("%s/%s/%s/%s_sce_nz_postQc%s_%sPca.png",
		 projDir, outDirBit, normPlotDirBit, setName, setSuf, typeNorm)
ggsave(filename=tmpFn, plot=p, type="cairo-png")
```


```r
tmpFn <- sprintf("%s/%s/%s_sce_nz_postQc%s_%sPca.png",
                 dirRel, normPlotDirBit, setName, setSuf, typeNorm)
knitr::include_graphics(tmpFn, auto_pdf = TRUE)
```

<img src="../Plots/Norm/GSM3872434_sce_nz_postQc_GSM3872434_logRawPca.png" width="1050" />

```r
rm(tmpFn)
```

<!--
Cell-wise RLE for the logRaw counts in the GSM3872434 set. Each cell is represented by a box plot showing the inter-quartile range in grey, wiskers colour-coded by Sample.Name and the median as a black circle. 
-->





### log CPM


```r
typeNorm <- "logCpm"

assay(sce, "logCpm") <- log2(calculateCPM(sce, size_factors = NULL) + 1)

logCpmPca <- runPCA(
  sce[,],
  exprs_values = "logCpm"
)
```

PCA plot for the logCpm counts in the GSM3872434 set.


```r
p <- plotPCA(
    logCpmPca,
    colour_by = "Sample.Name",
    size_by = "sum",
    shape_by = "source_name"
) + ggtitle(sprintf("PCA plot: %s", typeNorm))

# write plot to file:
tmpFn <- sprintf("%s/%s/%s/%s_sce_nz_postQc%s_%sPca.png",
		 projDir, outDirBit, normPlotDirBit, setName, setSuf, typeNorm)
ggsave(filename=tmpFn, plot=p, type="cairo-png")
```


```r
tmpFn <- sprintf("%s/%s/%s_sce_nz_postQc%s_%sPca.png",
                 dirRel, normPlotDirBit, setName, setSuf, typeNorm)
knitr::include_graphics(tmpFn, auto_pdf = TRUE)
```

<img src="../Plots/Norm/GSM3872434_sce_nz_postQc_GSM3872434_logCpmPca.png" width="1050" />

```r
rm(tmpFn)
```

<!--
Cell-wise RLE for the logCpm counts in the GSM3872434 set.
-->





### scran

Normalised counts are stored in the 'logcounts' assay


```r
typeNorm <- "scran"

# assay(sce, "logcounts")

scranPca <- runPCA(
  sceDeconv[,],
  exprs_values = "logcounts"
)
```

PCA plot for the 'scran' counts in the GSM3872434 set.


```r
p <- plotPCA(
    scranPca,
    colour_by = "Sample.Name",
    size_by = "sum",
    shape_by = "source_name"
) + ggtitle(sprintf("PCA plot: %s", typeNorm))

# write plot to file:
tmpFn <- sprintf("%s/%s/%s/%s_sce_nz_postQc%s_%sPca.png",
		 projDir, outDirBit, normPlotDirBit, setName, setSuf, typeNorm)
ggsave(filename=tmpFn, plot=p, type="cairo-png")
```


```r
tmpFn <- sprintf("%s/%s/%s_sce_nz_postQc%s_%sPca.png",
                 dirRel, normPlotDirBit, setName, setSuf, typeNorm)
knitr::include_graphics(tmpFn, auto_pdf = TRUE)
```

<img src="../Plots/Norm/GSM3872434_sce_nz_postQc_GSM3872434_scranPca.png" width="1050" />

```r
rm(tmpFn)
```

TSNE plot for the 'scran' counts in the GSM3872434 set.


```r
typeNorm <- "scran"

reducedDim(sceDeconv, "TSNE_scran") <- reducedDim(
  runTSNE(sceDeconv, exprs_values = "logcounts"),
  "TSNE"
)
```


```r
p <- plotReducedDim(
  sceDeconv,
  dimred = "TSNE_scran",
  colour_by = "Sample.Name",
  size_by = "sum",
  shape_by = "source_name"
) + ggtitle(sprintf("TSNE plot: %s", typeNorm))

# write plot to file:
tmpFn <- sprintf("%s/%s/%s/%s_sce_nz_postQc%s_%sTsne.png",
		 projDir, outDirBit, normPlotDirBit, setName, setSuf, typeNorm)
ggsave(filename=tmpFn, plot=p, type="cairo-png")
```


```r
tmpFn <- sprintf("%s/%s/%s_sce_nz_postQc%s_%sTsne.png",
                 dirRel, normPlotDirBit, setName, setSuf, typeNorm)
knitr::include_graphics(tmpFn, auto_pdf = TRUE)
```

<img src="../Plots/Norm/GSM3872434_sce_nz_postQc_GSM3872434_scranTsne.png" width="1050" />

```r
rm(tmpFn)
```

<!--
UMAP plot for the 'scran' counts in the GSM3872434 set.
-->







<!--
Cell-wise RLE for the 'scran' counts in the GSM3872434 set.
-->





### SCTransform


```r
typeNorm <- "sctrans"

reducedDim(sce, "PCA_sctrans_norm") <- reducedDim(
  runPCA(sce, exprs_values = "sctrans_norm"),
  "PCA"
)
```

PCA plot for the 'sctrans' counts in the GSM3872434 set.


```r
p <- plotReducedDim(
  sce,
  dimred = "PCA_sctrans_norm",
  colour_by = "Sample.Name",
  size_by = "sum",
  shape_by = "source_name"
) + ggtitle(sprintf("PCA plot: %s", typeNorm))

# write plot to file:
tmpFn <- sprintf("%s/%s/%s/%s_sce_nz_postQc%s_%sPca.png",
		 projDir, outDirBit, normPlotDirBit, setName, setSuf, typeNorm)
ggsave(filename=tmpFn, plot=p, type="cairo-png")
```


```r
tmpFn <- sprintf("%s/%s/%s_sce_nz_postQc%s_%sPca.png",
                 dirRel, normPlotDirBit, setName, setSuf, typeNorm)
knitr::include_graphics(tmpFn, auto_pdf = TRUE)
```

<img src="../Plots/Norm/GSM3872434_sce_nz_postQc_GSM3872434_sctransPca.png" width="1050" />

```r
rm(tmpFn)
```

TSNE plot for the sctrans counts in the GSM3872434 set.


```r
typeNorm <- "sctrans"

reducedDim(sce, "TSNE_sctrans_norm") <- reducedDim(
  runTSNE(sce, exprs_values = "sctrans_norm"),
  "TSNE"
)
```


```r
p <- plotReducedDim(
  sce,
  dimred = "TSNE_sctrans_norm",
  colour_by = "Sample.Name",
  size_by = "sum",
  shape_by = "source_name"
) + ggtitle(sprintf("TSNE plot: %s", typeNorm))

# write plot to file:
tmpFn <- sprintf("%s/%s/%s/%s_sce_nz_postQc%s_%sTsne.png",
		 projDir, outDirBit, normPlotDirBit, setName, setSuf, typeNorm)
ggsave(filename=tmpFn, plot=p, type="cairo-png")
```


```r
tmpFn <- sprintf("%s/%s/%s_sce_nz_postQc%s_%sTsne.png",
                 dirRel, normPlotDirBit, setName, setSuf, typeNorm)
knitr::include_graphics(tmpFn, auto_pdf = TRUE)
```

<img src="../Plots/Norm/GSM3872434_sce_nz_postQc_GSM3872434_sctransTsne.png" width="1050" />

```r
rm(tmpFn)
```

<!--
UMAP plot for the 'sctrans' counts in the GSM3872434 set.
-->







Cell-wise RLE for the sctrans counts in the GSM3872434 set.


```r
p <- plotRLE(
    sce,
    exprs_values = "sctrans_norm",
    colour_by = "Sample.Name"
) + ggtitle(sprintf("RLE plot: %s", typeNorm))

# write plot to file:
tmpFn <- sprintf("%s/%s/%s/%s_sce_nz_postQc%s_%sRle.png",
		 projDir, outDirBit, normPlotDirBit, setName, setSuf, typeNorm)
print("DEV"); print(getwd()); print(tmpFn)
```

```
## [1] "DEV"
```

```
## [1] "/ssd/personal/baller01/20200511_FernandesM_ME_crukBiSs2020/AnaWiSce/AnaCourse1/ScriptsCourse1"
```

```
## [1] "/ssd/personal/baller01/20200511_FernandesM_ME_crukBiSs2020/AnaWiSce/AnaCourse1/Plots/Norm/GSM3872434_sce_nz_postQc_GSM3872434_sctransRle.png"
```

```r
ggsave(filename=tmpFn, plot=p, type="cairo-png")
```


```r
tmpFn <- sprintf("%s/%s/%s_sce_nz_postQc%s_%sRle.png",
                 dirRel, normPlotDirBit, setName, setSuf, typeNorm)
knitr::include_graphics(tmpFn, auto_pdf = TRUE)
```

<img src="../Plots/Norm/GSM3872434_sce_nz_postQc_GSM3872434_sctransRle.png" width="1050" />

```r
rm(tmpFn)
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
##  [1] Cairo_1.5-12.2              BiocSingular_1.6.0         
##  [3] dplyr_1.0.6                 scran_1.18.7               
##  [5] scater_1.18.6               ggplot2_3.3.3              
##  [7] SingleCellExperiment_1.12.0 SummarizedExperiment_1.20.0
##  [9] Biobase_2.50.0              GenomicRanges_1.42.0       
## [11] GenomeInfoDb_1.26.7         IRanges_2.24.1             
## [13] S4Vectors_0.28.1            BiocGenerics_0.36.1        
## [15] MatrixGenerics_1.2.1        matrixStats_0.58.0         
## [17] knitr_1.33                 
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-7              sctransform_0.3.2.9005   
##  [3] tools_4.0.3               bslib_0.2.5              
##  [5] utf8_1.2.1                R6_2.5.0                 
##  [7] irlba_2.3.3               vipor_0.4.5              
##  [9] DBI_1.1.1                 colorspace_2.0-1         
## [11] withr_2.4.2               tidyselect_1.1.1         
## [13] gridExtra_2.3             compiler_4.0.3           
## [15] BiocNeighbors_1.8.2       isoband_0.2.4            
## [17] DelayedArray_0.16.3       labeling_0.4.2           
## [19] bookdown_0.22             sass_0.4.0               
## [21] scales_1.1.1              stringr_1.4.0            
## [23] digest_0.6.27             rmarkdown_2.8            
## [25] XVector_0.30.0            pkgconfig_2.0.3          
## [27] htmltools_0.5.1.1         parallelly_1.25.0        
## [29] sparseMatrixStats_1.2.1   limma_3.46.0             
## [31] highr_0.9                 rlang_0.4.11             
## [33] DelayedMatrixStats_1.12.3 jquerylib_0.1.4          
## [35] generics_0.1.0            farver_2.1.0             
## [37] jsonlite_1.7.2            BiocParallel_1.24.1      
## [39] RCurl_1.98-1.3            magrittr_2.0.1           
## [41] GenomeInfoDbData_1.2.4    scuttle_1.0.4            
## [43] Matrix_1.3-3              Rcpp_1.0.6               
## [45] ggbeeswarm_0.6.0          munsell_0.5.0            
## [47] fansi_0.4.2               viridis_0.6.1            
## [49] lifecycle_1.0.0           stringi_1.6.1            
## [51] yaml_2.2.1                edgeR_3.32.1             
## [53] MASS_7.3-54               zlibbioc_1.36.0          
## [55] Rtsne_0.15                plyr_1.8.6               
## [57] grid_4.0.3                listenv_0.8.0            
## [59] dqrng_0.3.0               crayon_1.4.1             
## [61] lattice_0.20-44           cowplot_1.1.1            
## [63] beachmat_2.6.4            locfit_1.5-9.4           
## [65] pillar_1.6.1              igraph_1.2.6             
## [67] future.apply_1.7.0        reshape2_1.4.4           
## [69] codetools_0.2-18          glue_1.4.2               
## [71] evaluate_0.14             vctrs_0.3.8              
## [73] png_0.1-7                 gtable_0.3.0             
## [75] purrr_0.3.4               future_1.21.0            
## [77] assertthat_0.2.1          xfun_0.23                
## [79] rsvd_1.0.5                viridisLite_0.4.0        
## [81] tibble_3.1.2              beeswarm_0.3.1           
## [83] globals_0.14.0            bluster_1.0.0            
## [85] statmod_1.4.36            ellipsis_0.3.2
```
