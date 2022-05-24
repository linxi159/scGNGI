# scGNGI
Missing value imputation with low-rank matrix completion in single-cell RNA-seq data by considering cell heterogeneity
![](https://github.com/linxi159/scGNGI/blob/master/figures/Figure_1_Final.jpg) 

## Description of each directory
data: the preprocessed data from simulated and real scRNA-seq data in GEO.

src: the implementation of missing value imputation with low-rank matrix completion

R: the utility of data preprocessing, comparison with different methods, and plot of experimental result.

results: the preprocessed results.

figures: all plots for scGNGI.


## How to setup

* Python (3.6 or later)

* numpy

* sklearn

* NVIDIA GPU + CUDA 11.50 + CuDNN v7.1

* scipy


## Quick example to use SRAFL
```
* train and test the model:

python ./src/main.py
```



