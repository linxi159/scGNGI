# scGNGI
Sequential reinforcement active feature learning for gene signature identification in renal cell carcinoma
![](https://github.com/linxi159/scGNGI/figures/Figure_1_Final.jpg) 

## Description of each directory
data: the preprocessed data from RCC in TCGA.

src: the implementation of Sequential reinforcement active feature learning.

R: the utility of data preprocessing, gene filtering, and plot of experimental result.

results: the selected gene signature and different cases.

figures: all plots for SRAFL.


## How to setup

* Python (3.6 or later)

* numpy

* sklearn

* NVIDIA GPU + CUDA 11.50 + CuDNN v7.1

* torch 0.8 or later

* collections

## Quick example to use SRAFL
```
* train and test the model:

python ./src/main.py
```



