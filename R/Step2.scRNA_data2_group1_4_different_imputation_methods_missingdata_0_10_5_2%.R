
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

i <- 10 # 10,5,2
dir.create(paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_scImpute"))
# dir.create(paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missing",i, "_scImpute"))
library(scImpute)
system.time(scimpute(count_path = paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing.csv"),
                     infile = "csv",           # format of input file
                     outfile = "csv",          # format of output file
                     out_dir = paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_scImpute"),  # full path to output directory
                     drop_thre = 0.5,          # threshold set on dropout probability
                     Kcluster = 4,
                     ncores = 5))
# system.time(scimpute(count_path = paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missing",i,".csv"),
#                      infile = "csv",           # format of input file
#                      outfile = "csv",          # format of output file
#                      out_dir = paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missing",i, "_scImpute"),  # full path to output directory
#                      drop_thre = 0.5,          # threshold set on dropout probability
#                      Kcluster = 4,
#                      ncores = 5))
# 

###############################################  
########### [2] running SAVER #################
############################################### 
setwd("~/r_workplace/8.genes_imputation/data")

i <- 10 # 10,5,2
dir.create(paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missing",i, "_SAVER"))
gene.expression <- read.csv(paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missing",i, ".csv"))
# dir.create(paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missing",i, "_SAVER"))
# gene.expression <- read.csv(paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missing",i, ".csv"))
rownames(gene.expression) <- gene.expression[,1]
gene.expression <- gene.expression[,-1]
library(SAVER)
gene.expression.missing <- gene.expression
system.time(saver_imp <- saver(gene.expression.missing))
save(saver_imp, file = paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_SAVER/scRNA_data2_group1_4_proprecessing_SAVER"))
#save(saver_imp, file = paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missing",i, "_SAVER/scRNA_data2_group1_4_proprecessing_missing",i,"_SAVER"))

saver_imp_estimate = saver_imp[["estimate"]]
write.csv(saver_imp_estimate, file = paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_SAVER/scRNA_data2_group1_4_proprecessing_SAVER_estimate.csv"))
#write.csv(saver_imp_estimate, file = paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missing",i,"_SAVER/scRNA_data2_group1_4_proprecessing_missing", i,"_SAVER_estimate.csv"))
saver_imp_se = saver_imp[["se"]]
write.csv(saver_imp_se, file = paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_SAVER/scRNA_data2_group1_4_proprecessing_SAVER_se.csv"))
#write.csv(saver_imp_se, file = paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missing",i,"_SAVER/scRNA_data2_group1_4_proprecessing_missing", i,"_SAVER_se.csv"))


###############################################
########### [3] running MAGIC #################
###############################################
setwd("~/r_workplace/8.genes_imputation/data")

i <- 10 #10,5,2
dir.create(paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_MAGIC"))
gene.expression <- read.csv(paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing.csv"))
# dir.create(paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missing",i, "_MAGIC"))
# gene.expression <- read.csv(paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missing",i, ".csv"))
rownames(gene.expression) <- gene.expression[,1]
gene.expression <- gene.expression[,-1]
gene.expression.t <- t(gene.expression)
gene.expression.t <- as.data.frame(gene.expression.t)

library(Rmagic)
library(ggplot2)

data_MAGIC <- magic(gene.expression.t, genes="all_genes")
data_MAGIC_t <- t(data_MAGIC[["result"]])
data_MAGIC_t <- as.data.frame(data_MAGIC_t)

write.csv(data_MAGIC_t, file = paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_MAGIC/scRNA_data2_group1_4_proprecessing_MAGIC.csv"))
#write.csv(data_MAGIC_t, file = paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missing",i,"_MAGIC/scRNA_data2_group1_4_proprecessing_missing", i,"_MAGIC.csv"))


####################################################
################ [4] running DrImpute) #############
####################################################
# DrImpute
setwd("~/r_workplace/8.genes_imputation/data")
library(ADImpute)
i <- 10 #10,5,2
dir.create(paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_DrImpute"))
gene.expression <- read.csv(paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing.csv"))
# dir.create(paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missing",i, "_DrImpute"))
# gene.expression <- read.csv(paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missing",i, ".csv"))
rownames(gene.expression) <- gene.expression[,1]
gene.expression <- gene.expression[,-1]
data_ <- as.matrix(gene.expression)

imputed <- Impute(data = data_ , do = c("DrImpute"), cores = 2)
imputed <- as.data.frame(imputed[["DrImpute"]]) 
write.csv(imputed, file = paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_DrImpute/scRNA_data2_group1_4_proprecessing_DrImpute.csv"))
#write.csv(imputed, file = paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missing",i,"_DrImpute/scRNA_data2_group1_4_proprecessing_missing", i,"_DrImpute.csv"))


################################################
############ [5] running ALRA ##################
################################################
setwd("~/r_workplace/8.genes_imputation/data")
source('alra.R')
i <- 10 #10,5,2
dir.create(paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_ALRA"))
gene.expression <- read.csv(paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing.csv"))
# dir.create(paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missing",i, "_ALRA"))
# gene.expression <- read.csv(paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missing",i, ".csv"))
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
write.csv(imputed, file = paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_ALRA/scRNA_data2_group1_4_proprecessing_ALRA.csv"))
#write.csv(imputed, file = paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missing",i,"_ALRA/scRNA_data2_group1_4_proprecessing_missing", i,"_ALRA.csv"))


###############################################
########## [6] running DeepImpute #############
###############################################
setwd("~/r_workplace/8.genes_imputation/data")

i <- 10 # 10,5,2)
dir.create(paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missing",i, "_DeepImput"))

# python 执行
# import deepimpute
# from deepimpute.multinet import MultiNet
# import pandas as pd
#
# arg_train = True # False
# df = pd.read_csv('./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missingi.csv', index_col=0) # dimension = (genes x cells)
# data = pd.DataFrame(df.values.T, index=df.columns, columns=df.index) # dimension = (cells x genes)
# model = MultiNet()
# if arg_train == True:
#   model.fit(data)
# else:
#   model.load()
#
# imputed = model.predict(data)
# imputed_ = pd.DataFrame(imputed.values.T, index=imputed.columns, columns=imputed.index) # dimension = (genes x cells)
# # save
# imputed_.to_csv('./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missingi_DeepImpute/scRNA_data2_group1_4_proprecessing_missingi_DeepImpute.csv')


###############################################
##### [7] running scGNGI(proposedMethods) #####
###############################################

setwd("~/r_workplace/8.genes_imputation/data")

i <- 10 #10,5,2
dir.create(paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missing",i, "_scGNGI"))

# python 执行



