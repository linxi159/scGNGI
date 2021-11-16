
###############################################
########### [1] running scImpute ##############
###############################################
setwd("~/r_workplace/8.genes_imputation/data")

## WORKAROUND: https://github.com/rstudio/rstudio/issues/6692
## Revert to 'sequential' setup of PSOCK cluster in RStudio Console on macOS and R 4.0.0
if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) &&
    Sys.info()["sysname"] == "Darwin" && getRversion() >= "4.0.0") {
  parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
}

i <- 35 # 10,35
dir.create(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i, "_scImpute"))
library(scImpute)
system.time(scimpute(count_path = paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,".csv"),
                     infile = "csv",           # format of input file
                     outfile = "csv",          # format of output file
                     out_dir = paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i, "_scImpute"),  # full path to output directory
                     drop_thre = 0.5,          # threshold set on dropout probability
                     Kcluster = 3,
                     ncores = 5))


###############################################  
########### [2] running SAVER #################
############################################### 
rm(list = ls())
setwd("~/r_workplace/8.genes_imputation/data")

i <- 35 # 10,35
dir.create(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i, "_SAVER"))
gene.expression <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i, ".csv"))
rownames(gene.expression) <- gene.expression[,1]
gene.expression <- gene.expression[,-1]
library(SAVER)
gene.expression.missing <- gene.expression
system.time(saver_imp <- saver(gene.expression.missing))
save(saver_imp, file = paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i, "_SAVER/scRNA_data3_group1_3_proprecessing_missing",i,"_SAVER"))

saver_imp_estimate = saver_imp[["estimate"]]
write.csv(saver_imp_estimate, file = paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,"_SAVER/scRNA_data3_group1_3_proprecessing_missing", i,"_SAVER_estimate.csv"))
saver_imp_se = saver_imp[["se"]]
write.csv(saver_imp_se, file = paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,"_SAVER/scRNA_data3_group1_3_proprecessing_missing", i,"_SAVER_se.csv"))


###############################################
########### [3] running MAGIC #################
###############################################
setwd("~/r_workplace/8.genes_imputation/data")

i <- 35 # 10,35
dir.create(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i, "_MAGIC"))
gene.expression <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i, ".csv"))
rownames(gene.expression) <- gene.expression[,1]
gene.expression <- gene.expression[,-1]
gene.expression.t <- t(gene.expression)
gene.expression.t <- as.data.frame(gene.expression.t)

library(Rmagic)
library(ggplot2)

data_MAGIC <- magic(gene.expression.t, genes="all_genes")
data_MAGIC_t <- t(data_MAGIC[["result"]])
data_MAGIC_t <- as.data.frame(data_MAGIC_t)

write.csv(data_MAGIC_t, file = paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,"_MAGIC/scRNA_data3_group1_3_proprecessing_missing", i,"_MAGIC.csv"))


####################################################
################ [4] running DrImpute) #############
####################################################
# DrImpute
setwd("~/r_workplace/8.genes_imputation/data")
library(ADImpute)
i <- 35 # 10,35
dir.create(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i, "_DrImpute"))
gene.expression <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i, ".csv"))
rownames(gene.expression) <- gene.expression[,1]
gene.expression <- gene.expression[,-1]
data_ <- as.matrix(gene.expression)

imputed <- Impute(data = data_ , do = c("DrImpute"), cores = 2)
imputed <- as.data.frame(imputed[["DrImpute"]]) 
write.csv(imputed, file = paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,"_DrImpute/scRNA_data3_group1_3_proprecessing_missing", i,"_DrImpute.csv"))


################################################
############ [5] running ALRA ##################
################################################
setwd("~/r_workplace/8.genes_imputation/data")
source('alra.R')
i <- 35 # 10,15,20,25,30,35
dir.create(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i, "_ALRA"))
gene.expression <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i, ".csv"))
rownames(gene.expression) <- gene.expression[,1]
gene.expression <- gene.expression[,-1]
gene.expression.t <- t(gene.expression)
gene.expression.t <- log2(gene.expression.t+1)

# Let A_norm be a normalized expression matrix where cells are rows and genes are columns.
# We use library and log normalization, but other approaches may also work well.
result.completed <- alra(gene.expression.t)
imputed <- result.completed[[3]]

imputed <- (2^(imputed)-1)
imputed <- t(imputed)
colnames(imputed) <- rownames(gene.expression.t)
imputed <- as.data.frame(imputed)
write.csv(imputed, file = paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,"_ALRA/scRNA_data3_group1_3_proprecessing_missing", i,"_ALRA.csv"))


###############################################
##### [6] running scGNGI(proposedMethods) #####
###############################################

setwd("~/r_workplace/8.genes_imputation/data")

i <- 10 # 10,15,20,25,30,35
dir.create(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i, "_scGNGI"))

# python 执行



