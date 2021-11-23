## Inferring Gene Regulatory Network (GRN) using partial information decomposition and context (PIDC)

### 1. Installing PIDC from github

The PIDC package is implemented with R and deposited at Github: https://github.com/hmutpw/PIDC. To using the PIDC, you must install the following dependent packages: **purrr**, **parallel**, **pbapply**, **reshape2**, **minet**.

```{r}
install.packages("devtools")
install.packages("purrr")
install.packages("parallel")
install.packages("pbapply")
install.packages("reshape2")
install.packages("minet")

devtools::install_github("hmutpw/PIDC")
```

### 2. Running PIDC with test data 

The input matrix of PIDC should be non-negative values such as raw counts, UMIs, TPMs or FPKMs. Log-scaled expressed value is suggested. Each row of the matrix represent the feature such as gene or isoform, while each column represent cell. Here we use test data with 50 genes and 4488 cells in **PIDC** package.

```{r}
library(PIDC)
data(expMat)
# 1. filtering genes with expression (UMIs >= 1) in at least 10 cells
expMat <- expMat[apply(expMat,1,function(x){length(x[which(x>=1)])>=10}),]

# 2. Runing PIDC using UMIs without log scaled.
PIDC_grn <- PIDC(expMat = expMat, logScale = TRUE, ncores = 1)
```

The output of PIDC is a square matrix with weighted between each gene pairs.

```{r}
PIDC_grn[1:5,1:5]
```

**Note**:The running time heavily depend on the number of genes in your input matrix, for genes over 5,000, we strongly suggest you to use multi-core computers. Generally, we recommend users **NOT** to set regulators and targets to get accurate results, because the PIDC algorithm calculated weights based on triple gene pairs, absence of regulators or targets will decrease accuracy.

### 3. Generating gene regulons from GRNs

To remove weakest edges from gene regulatory network and keep the strong connected edges for further analysis. We using `aracne` to generate the gene regulons from GRNs. The input of `matToNet` is a square matrix with positive values.

```{r}
PIDC_net <- matToNet(weightMat = PIDC_grn, methods = "aracne")
```

The final output is a three-column network with weights:

```{r}
head(PIDC_net)
```

### 4. sessionInfo

```{r}
sessionInfo()
```

