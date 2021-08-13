# PIDC
Inferring Gene Regulatory Network (GRN) using partial information decomposition and context (PIDC) method.
#### 1. Install PIDC

```r
install.packages("devtools")
library(devtools)
install_github("hmutpw/PIDC")
library(PIDC)
```

#### 2. Running PIDC with test data

The input matrix of PIDC should be non-negative raw count, UMI, TPM or log-scaled expressed values. Rows represent for genes and columns for samples (or cells).

```R
data(expMat)
PIDC_res <- PIDC(expMat)
```

The running time heavily depend on the number of genes in expression matrix, for genes over 5,000,  we strongly suggest you to use multi-core computers. Meanwhile, we recommend users **NOT** to set  *regulators* and *targets* to get better results, because the PIDC will calculate weights based on  triple gene pairs.

#### 3. Generating regulatory networks

```R
PIDC_net <- matToNet(mat=PIDC_res)
```

Here, we use different methods to remove weakest edges and keep the strong connected edges for further analysis, for example, transcription factor regulon activity analysis usinng **[scATFR](https://github.com/hmutpw/scATFR)**.  



