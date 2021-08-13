# Inferring Gene Regulatory Network (GRN) using partial information decomposition and context (PIDC) method

#### 1. Install PIDC

```r
install.packages("devtools")
library(devtools)
install_github("hmutpw/PIDC")
library(PIDC)
```

#### 2. Running PIDC with test data

```R
data(expMat)
PIDC_res <- PIDC(expMat)
```

The running time heavily depend on the number of regulators and targets, for genes over 5000,  we strongly suggest you to use multi-cores computers. Meanwhile, we suggest you do not set the regulators and targets which may perform better since the PIDC calculating the regulatory relationships between every triple gene pairs.

#### 3. Generating regulatory networks

```R
PIDC_net <- matToNet(mat=PIDC_res)
```

Here, we use different methods to remove weakest edges and keep the strong connected edges for further analysis, for example, transcription factor regulon activity analysis usinng **[scATFR](https://github.com/hmutpw/scATFR)**.  
