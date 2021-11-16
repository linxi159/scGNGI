
###############################################
########### Visualization  Examples ###########
###############################################

################################################# i = 10  ###################################################
#################################################### [0] ####################################################
# load data
rm(list=ls())
setwd("~/r_workplace/8.genes_imputation/data")

#################################################### [1] ####################################################
i = 10 
# raw_data [1]
df_raw <- read.csv("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing.csv") 

rownames(df_raw) <- df_raw[,1]
df_raw <- df_raw[,-1]

method_raw = "Groud Truth" #[1.0][2.0]

df <- df_raw[,] 

df <- t(df)
df <- as.data.frame(df)
class <- substring(rownames(df), 1,6) # scRNA_data3_group1_3
label <- vector()
for (j in c(1:length(class)) ) {
  label[j] <- class[[j]][1]
}
df$label <- label

# 加载tsne包
library(Rtsne)
library(umap)
colors = rainbow(length(unique(df$label)))
names(colors) = unique(df$label)

# 使用tsne或者Rtsne函数进行tSNE降维分析
set.seed(42)
tsne_unique <- unique(df) # Remove duplicates
tsne_matrix <- as.matrix(tsne_unique[,1:(length(df)-1)])

Sys.time()
tsne_out <- Rtsne(tsne_matrix,dims=2,perplexity=30,theta=0.0,pca = FALSE) # Run TSNE
tsne_raw <- tsne_out$Y

umap_out <- umap(tsne_matrix)
umap_raw <- umap_out$layout
Sys.time()

#################################################### [2] ####################################################
# raw_data_missing10 [2]
df_raw <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,".csv")) 

rownames(df_raw) <- df_raw[,1]
df_raw <- df_raw[,-1]

method_raw_missing = "raw_missing10" 

df <- df_raw[,] 

df <- t(df)
df <- as.data.frame(df)
class <- substring(rownames(df), 1,6) # scRNA_data1_group1_5 // scRNA_data2_group1_4
label <- vector()
for (j in c(1:length(class)) ) {
  label[j] <- class[[j]][1]
}
df$label <- label

# 加载tsne包
library(Rtsne)
library(umap)
colors = rainbow(length(unique(df$label)))
names(colors) = unique(df$label)

# 使用tsne或者Rtsne函数进行tSNE降维分析
set.seed(42)
tsne_unique <- unique(df) # Remove duplicates
tsne_matrix <- as.matrix(tsne_unique[,1:(length(df)-1)])

Sys.time()
tsne_out <- Rtsne(tsne_matrix,dims=2,perplexity=30,theta=0.0,pca = FALSE) # Run TSNE
tsne_raw_missing <- tsne_out$Y

umap_out <- umap(tsne_matrix)
umap_raw_missing <- umap_out$layout
Sys.time()

#################################################### [3] ####################################################
method = "MAGIC"
method_MAGIC = "MAGIC"
file_name = "MAGIC.csv"

# missing_impute_data_with_methods [2.0]
df_method <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,"_",method, "/scRNA_data3_group1_3_proprecessing_missing",i,"_",file_name))

rownames(df_method) <- df_method[,1]
df_method <- df_method[,-1]

# mask missing_zero_space 1
df_mask <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,"_mask.csv"))
rownames(df_mask) <- df_mask[,1]
df_mask <- df_mask[,-1]
df_mask <- as.matrix(df_mask)

df <- df_raw[,] + df_method*df_mask

df <- t(df)
df <- as.data.frame(df)
class <- substring(rownames(df), 1,6) # scRNA_data1_group1_5 // scRNA_data2_group1_4
label <- vector()
for (j in c(1:length(class)) ) {
  label[j] <- class[[j]][1]
}
df$label <- label

# 加载tsne包
library(Rtsne)
library(umap)
colors = rainbow(length(unique(df$label)))
names(colors) = unique(df$label)

# 使用tsne或者Rtsne函数进行tSNE降维分析
set.seed(42)
tsne_unique <- unique(df) # Remove duplicates
tsne_matrix <- as.matrix(tsne_unique[,1:(length(df)-1)])

Sys.time()
tsne_out <- Rtsne(tsne_matrix,dims=2,perplexity=30,theta=0.0,pca = FALSE) # Run TSNE
tsne_MAGIC <- tsne_out$Y

umap_out <- umap(tsne_matrix)
umap_MAGIC <- umap_out$layout
Sys.time()

#################################################### [4] ####################################################
method = "SAVER"
method_SAVER = "SAVER"
file_name = "SAVER_estimate.csv"

# missing_impute_data_with_methods [2.0]
df_method <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,"_",method, "/scRNA_data3_group1_3_proprecessing_missing",i,"_",file_name))

rownames(df_method) <- df_method[,1]
df_method <- df_method[,-1]

# mask missing_zero_space 1
df_mask <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,"_mask.csv"))
rownames(df_mask) <- df_mask[,1]
df_mask <- df_mask[,-1]
df_mask <- as.matrix(df_mask)

df <- df_raw[,] + df_method*df_mask

df <- t(df)
df <- as.data.frame(df)
class <- substring(rownames(df), 1,6) # scRNA_data1_group1_5 // scRNA_data2_group1_4
label <- vector()
for (j in c(1:length(class)) ) {
  label[j] <- class[[j]][1]
}
df$label <- label

# 加载tsne包
library(Rtsne)
library(umap)
colors = rainbow(length(unique(df$label)))
names(colors) = unique(df$label)

# 使用tsne或者Rtsne函数进行tSNE降维分析
set.seed(42)
tsne_unique <- unique(df) # Remove duplicates
tsne_matrix <- as.matrix(tsne_unique[,1:(length(df)-1)])

Sys.time()
tsne_out <- Rtsne(tsne_matrix,dims=2,perplexity=30,theta=0.0,pca = FALSE) # Run TSNE
tsne_SAVER <- tsne_out$Y

umap_out <- umap(tsne_matrix)
umap_SAVER <- umap_out$layout
Sys.time()

#################################################### [5] ####################################################
method = "DrImpute"
method_DrImpute = "DrImpute"
file_name = "DrImpute.csv"

# missing_impute_data_with_methods [2.0]
df_method <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,"_",method, "/scRNA_data3_group1_3_proprecessing_missing",i,"_",file_name))

rownames(df_method) <- df_method[,1]
df_method <- df_method[,-1]

# mask missing_zero_space 1
df_mask <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,"_mask.csv"))
rownames(df_mask) <- df_mask[,1]
df_mask <- df_mask[,-1]
df_mask <- as.matrix(df_mask)

df <- df_raw[,] + df_method*df_mask

df <- t(df)
df <- as.data.frame(df)
class <- substring(rownames(df), 1,6) # scRNA_data1_group1_5 // scRNA_data2_group1_4
label <- vector()
for (j in c(1:length(class)) ) {
  label[j] <- class[[j]][1]
}
df$label <- label

# 加载tsne包
library(Rtsne)
library(umap)
colors = rainbow(length(unique(df$label)))
names(colors) = unique(df$label)

# 使用tsne或者Rtsne函数进行tSNE降维分析
set.seed(42)
tsne_unique <- unique(df) # Remove duplicates
tsne_matrix <- as.matrix(tsne_unique[,1:(length(df)-1)])

Sys.time()
tsne_out <- Rtsne(tsne_matrix,dims=2,perplexity=30,theta=0.0,pca = FALSE) # Run TSNE
tsne_DrImpute <- tsne_out$Y

umap_out <- umap(tsne_matrix)
umap_DrImpute <- umap_out$layout
Sys.time()

#################################################### [6] ####################################################
method = "scGNGI"
method_scGNGI = "scGNGI"
file_name = "scGNGI_iter1_r3.csv"

# missing_impute_data_with_methods [2.0]
df_method <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,"_",method, "/scRNA_data3_group1_3_proprecessing_missing",i,"_",file_name))

rownames(df_method) <- df_method[,1]
df_method <- df_method[,-1]

# mask missing_zero_space 1
df_mask <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,"_mask.csv"))
rownames(df_mask) <- df_mask[,1]
df_mask <- df_mask[,-1]
df_mask <- as.matrix(df_mask)

df <- df_raw[,] + df_method*df_mask

df <- t(df)
df <- as.data.frame(df)
class <- substring(rownames(df), 1,6) # scRNA_data1_group1_5 // scRNA_data2_group1_4
label <- vector()
for (j in c(1:length(class)) ) {
  label[j] <- class[[j]][1]
}
df$label <- label

# 加载tsne包
library(Rtsne)
library(umap)
colors = rainbow(length(unique(df$label)))
names(colors) = unique(df$label)

# 使用tsne或者Rtsne函数进行tSNE降维分析
set.seed(42)
tsne_unique <- unique(df) # Remove duplicates
tsne_matrix <- as.matrix(tsne_unique[,1:(length(df)-1)])

Sys.time()
tsne_out <- Rtsne(tsne_matrix,dims=2,perplexity=30,theta=0.0,pca = FALSE) # Run TSNE
tsne_scGNGI <- tsne_out$Y

umap_out <- umap(tsne_matrix)
umap_scGNGI <- umap_out$layout
Sys.time()

save(tsne_raw, tsne_raw_missing, tsne_MAGIC, tsne_SAVER, tsne_DrImpute, tsne_scGNGI, umap_raw, umap_raw_missing, umap_MAGIC, umap_SAVER, umap_DrImpute, umap_scGNGI, file = paste0("./results/visualization/scRNA_data3_group1_3/tsne_umap_results.summary"))

#################################################### [7] ####################################################
setwd("~/r_workplace/8.genes_imputation/data")
load(paste0("./results/visualization/scRNA_data3_group1_3/tsne_umap_results.summary"))

#### plots tsne####
# 分割屏幕
split.screen(c(2, 3))

# 绘制第一个图
screen(1)
plot(tsne_raw_missing,col=colors[df$label],pch=16,
     xlab = "t-SNE1",ylab = "t-SNE2",main = method_raw_missing, cex.axis=0.8,cex.lab=0.8,cex.main=0.8)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 绘制第二个图
screen(2)
plot(tsne_MAGIC,col=colors[df$label],pch=16,
     xlab = "t-SNE1",ylab = "t-SNE2",main = method_MAGIC,cex.axis=0.8,cex.lab=0.8,cex.main=0.8)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 绘制第三个图
screen(3)
plot(tsne_SAVER,col=colors[df$label],pch=16,
     xlab = "t-SNE1",ylab = "t-SNE2",main = method_SAVER,cex.axis=0.8,cex.lab=0.8,cex.main=0.8)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 绘制第四个图
screen(4)
plot(tsne_DrImpute,col=colors[df$label],pch=16,
     xlab = "t-SNE1",ylab = "t-SNE2",main = method_DrImpute,cex.axis=0.8,cex.lab=0.8,cex.main=0.8)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 绘制第五个图
screen(5)
plot(tsne_scGNGI,col=colors[df$label],pch=16,
     xlab = "t-SNE1",ylab = "t-SNE2",main = method_scGNGI,cex.axis=0.8,cex.lab=0.8,cex.main=0.8)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 绘制第六个图
screen(6)
plot(tsne_raw,col=colors[df$label],pch=16,
     xlab = "t-SNE1",ylab = "t-SNE2",main = method_raw,cex.axis=0.8,cex.lab=0.8,cex.main=0.8)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 结束绘图
dev.off()

#### plots umap####
# 分割屏幕
split.screen(c(2, 3))

# 绘制第一个图
screen(1)
plot(umap_raw_missing,col=colors[df$label],pch=16,
     xlab = "UMAP_1",ylab = "UMAP_2",main = method_raw_missing,cex.axis=0.8,cex.lab=0.8,cex.main=0.8)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 绘制第二个图
screen(2)
plot(umap_MAGIC,col=colors[df$label],pch=16,
     xlab = "UMAP_1",ylab = "UMAP_2",main = method_MAGIC,cex.axis=0.8,cex.lab=0.8,cex.main=0.8)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 绘制第三个图
screen(3)
plot(umap_SAVER,col=colors[df$label],pch=16,
     xlab = "UMAP_1",ylab = "UMAP_2",main = method_SAVER,cex.axis=0.8,cex.lab=0.8,cex.main=0.8)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 绘制第四个图
screen(4)
plot(umap_DrImpute,col=colors[df$label],pch=16,
     xlab = "UMAP_1",ylab = "UMAP_2",main = method_DrImpute,cex.axis=0.8,cex.lab=0.8,cex.main=0.8)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 绘制第五个图
screen(5)
plot(umap_scGNGI,col=colors[df$label],pch=16,
     xlab = "UMAP_1",ylab = "UMAP_2",main = method_scGNGI,cex.axis=0.8,cex.lab=0.8,cex.main=0.8)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 绘制第六个图
screen(6)
plot(umap_raw,col=colors[df$label],pch=16,
     xlab = "UMAP_1",ylab = "UMAP_2",main = method_raw,cex.axis=0.8,cex.lab=0.8,cex.main=0.8)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 结束绘图
dev.off()

# 保存大小 width X height: pdf- 12 8,  eps- 1200 800


################################################# i = 35  ################################################
#################################################### [0] ####################################################
# load data
rm(list=ls())
setwd("~/r_workplace/8.genes_imputation/data")

#################################################### [1] ####################################################
i = 35 
# raw_data [1]
df_raw <- read.csv("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing.csv") 

rownames(df_raw) <- df_raw[,1]
df_raw <- df_raw[,-1]

method_raw = "Groud Truth" #[1.0][2.0]

df <- df_raw[,] 

df <- t(df)
df <- as.data.frame(df)
class <- substring(rownames(df), 1,6) # scRNA_data3_group1_3
label <- vector()
for (j in c(1:length(class)) ) {
  label[j] <- class[[j]][1]
}
df$label <- label

# 加载tsne包
library(Rtsne)
library(umap)
colors = rainbow(length(unique(df$label)))
names(colors) = unique(df$label)

# 使用tsne或者Rtsne函数进行tSNE降维分析
set.seed(42)
tsne_unique <- unique(df) # Remove duplicates
tsne_matrix <- as.matrix(tsne_unique[,1:(length(df)-1)])

Sys.time()
tsne_out <- Rtsne(tsne_matrix,dims=2,perplexity=30,theta=0.0,pca = FALSE) # Run TSNE
tsne_raw <- tsne_out$Y

umap_out <- umap(tsne_matrix)
umap_raw <- umap_out$layout
Sys.time()

#################################################### [2] ####################################################
# raw_data_missing10 [2]
df_raw <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,".csv")) 

rownames(df_raw) <- df_raw[,1]
df_raw <- df_raw[,-1]

method_raw_missing = "raw_missing35" 

df <- df_raw[,] 

df <- t(df)
df <- as.data.frame(df)
class <- substring(rownames(df), 1,6) # scRNA_data1_group1_5 // scRNA_data2_group1_4
label <- vector()
for (j in c(1:length(class)) ) {
  label[j] <- class[[j]][1]
}
df$label <- label

# 加载tsne包
library(Rtsne)
library(umap)
colors = rainbow(length(unique(df$label)))
names(colors) = unique(df$label)

# 使用tsne或者Rtsne函数进行tSNE降维分析
set.seed(42)
tsne_unique <- unique(df) # Remove duplicates
tsne_matrix <- as.matrix(tsne_unique[,1:(length(df)-1)])

Sys.time()
tsne_out <- Rtsne(tsne_matrix,dims=2,perplexity=30,theta=0.0,pca = FALSE) # Run TSNE
tsne_raw_missing <- tsne_out$Y

umap_out <- umap(tsne_matrix)
umap_raw_missing <- umap_out$layout
Sys.time()

#################################################### [3] ####################################################
method = "MAGIC"
method_MAGIC = "MAGIC"
file_name = "MAGIC.csv"

# missing_impute_data_with_methods [2.0]
df_method <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,"_",method, "/scRNA_data3_group1_3_proprecessing_missing",i,"_",file_name))

rownames(df_method) <- df_method[,1]
df_method <- df_method[,-1]

# mask missing_zero_space 1
df_mask <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,"_mask.csv"))
rownames(df_mask) <- df_mask[,1]
df_mask <- df_mask[,-1]
df_mask <- as.matrix(df_mask)

df <- df_raw[,] + df_method*df_mask

df <- t(df)
df <- as.data.frame(df)
class <- substring(rownames(df), 1,6) # scRNA_data1_group1_5 // scRNA_data2_group1_4
label <- vector()
for (j in c(1:length(class)) ) {
  label[j] <- class[[j]][1]
}
df$label <- label

# 加载tsne包
library(Rtsne)
library(umap)
colors = rainbow(length(unique(df$label)))
names(colors) = unique(df$label)

# 使用tsne或者Rtsne函数进行tSNE降维分析
set.seed(42)
tsne_unique <- unique(df) # Remove duplicates
tsne_matrix <- as.matrix(tsne_unique[,1:(length(df)-1)])

Sys.time()
tsne_out <- Rtsne(tsne_matrix,dims=2,perplexity=30,theta=0.0,pca = FALSE) # Run TSNE
tsne_MAGIC <- tsne_out$Y

umap_out <- umap(tsne_matrix)
umap_MAGIC <- umap_out$layout
Sys.time()

#################################################### [4] ####################################################
method = "scImpute"
method_scImpute = "scImpute"
file_name = "scImputescimpute_count.csv"

# missing_impute_data_with_methods [2.0]
df_method <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,"_",method, "/scRNA_data3_group1_3_proprecessing_missing",i,"_",file_name))

rownames(df_method) <- df_method[,1]
df_method <- df_method[,-1]

# mask missing_zero_space 1
df_mask <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,"_mask.csv"))
rownames(df_mask) <- df_mask[,1]
df_mask <- df_mask[,-1]
df_mask <- as.matrix(df_mask)

df <- df_raw[,] + df_method*df_mask

df <- t(df)
df <- as.data.frame(df)
class <- substring(rownames(df), 1,6) # scRNA_data1_group1_5 // scRNA_data2_group1_4
label <- vector()
for (j in c(1:length(class)) ) {
  label[j] <- class[[j]][1]
}
df$label <- label

# 加载tsne包
library(Rtsne)
library(umap)
colors = rainbow(length(unique(df$label)))
names(colors) = unique(df$label)

# 使用tsne或者Rtsne函数进行tSNE降维分析
set.seed(42)
tsne_unique <- unique(df) # Remove duplicates
tsne_matrix <- as.matrix(tsne_unique[,1:(length(df)-1)])

Sys.time()
tsne_out <- Rtsne(tsne_matrix,dims=2,perplexity=30,theta=0.0,pca = FALSE) # Run TSNE
tsne_scImpute <- tsne_out$Y

umap_out <- umap(tsne_matrix)
umap_scImpute <- umap_out$layout
Sys.time()

#################################################### [5] ####################################################
method = "SAVER"
method_SAVER = "SAVER"
file_name = "SAVER_estimate.csv"

# missing_impute_data_with_methods [2.0]
df_method <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,"_",method, "/scRNA_data3_group1_3_proprecessing_missing",i,"_",file_name))

rownames(df_method) <- df_method[,1]
df_method <- df_method[,-1]

# mask missing_zero_space 1
df_mask <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,"_mask.csv"))
rownames(df_mask) <- df_mask[,1]
df_mask <- df_mask[,-1]
df_mask <- as.matrix(df_mask)

df <- df_raw[,] + df_method*df_mask

df <- t(df)
df <- as.data.frame(df)
class <- substring(rownames(df), 1,6) # scRNA_data1_group1_5 // scRNA_data2_group1_4
label <- vector()
for (j in c(1:length(class)) ) {
  label[j] <- class[[j]][1]
}
df$label <- label

# 加载tsne包
library(Rtsne)
library(umap)
colors = rainbow(length(unique(df$label)))
names(colors) = unique(df$label)

# 使用tsne或者Rtsne函数进行tSNE降维分析
set.seed(42)
tsne_unique <- unique(df) # Remove duplicates
tsne_matrix <- as.matrix(tsne_unique[,1:(length(df)-1)])

Sys.time()
tsne_out <- Rtsne(tsne_matrix,dims=2,perplexity=30,theta=0.0,pca = FALSE) # Run TSNE
tsne_SAVER <- tsne_out$Y

umap_out <- umap(tsne_matrix)
umap_SAVER <- umap_out$layout
Sys.time()

#################################################### [6] ####################################################
method = "DrImpute"
method_DrImpute = "DrImpute"
file_name = "DrImpute.csv"

# missing_impute_data_with_methods [2.0]
df_method <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,"_",method, "/scRNA_data3_group1_3_proprecessing_missing",i,"_",file_name))

rownames(df_method) <- df_method[,1]
df_method <- df_method[,-1]

# mask missing_zero_space 1
df_mask <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,"_mask.csv"))
rownames(df_mask) <- df_mask[,1]
df_mask <- df_mask[,-1]
df_mask <- as.matrix(df_mask)

df <- df_raw[,] + df_method*df_mask

df <- t(df)
df <- as.data.frame(df)
class <- substring(rownames(df), 1,6) # scRNA_data1_group1_5 // scRNA_data2_group1_4
label <- vector()
for (j in c(1:length(class)) ) {
  label[j] <- class[[j]][1]
}
df$label <- label

# 加载tsne包
library(Rtsne)
library(umap)
colors = rainbow(length(unique(df$label)))
names(colors) = unique(df$label)

# 使用tsne或者Rtsne函数进行tSNE降维分析
set.seed(42)
tsne_unique <- unique(df) # Remove duplicates
tsne_matrix <- as.matrix(tsne_unique[,1:(length(df)-1)])

Sys.time()
tsne_out <- Rtsne(tsne_matrix,dims=2,perplexity=30,theta=0.0,pca = FALSE) # Run TSNE
tsne_DrImpute <- tsne_out$Y

umap_out <- umap(tsne_matrix)
umap_DrImpute <- umap_out$layout
Sys.time()

#################################################### [7] ####################################################
method = "ALRA"
method_ALRA = "ALRA"
file_name = "ALRA.csv"

# missing_impute_data_with_methods [2.0]
df_method <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,"_",method, "/scRNA_data3_group1_3_proprecessing_missing",i,"_",file_name))

rownames(df_method) <- df_method[,1]
df_method <- df_method[,-1]

# mask missing_zero_space 1
df_mask <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,"_mask.csv"))
rownames(df_mask) <- df_mask[,1]
df_mask <- df_mask[,-1]
df_mask <- as.matrix(df_mask)

df <- df_raw[,] + df_method*df_mask

df <- t(df)
df <- as.data.frame(df)
class <- substring(rownames(df), 1,6) # scRNA_data1_group1_5 // scRNA_data2_group1_4
label <- vector()
for (j in c(1:length(class)) ) {
  label[j] <- class[[j]][1]
}
df$label <- label

# 加载tsne包
library(Rtsne)
library(umap)
colors = rainbow(length(unique(df$label)))
names(colors) = unique(df$label)

# 使用tsne或者Rtsne函数进行tSNE降维分析
set.seed(42)
tsne_unique <- unique(df) # Remove duplicates
tsne_matrix <- as.matrix(tsne_unique[,1:(length(df)-1)])

Sys.time()
tsne_out <- Rtsne(tsne_matrix,dims=2,perplexity=30,theta=0.0,pca = FALSE) # Run TSNE
tsne_ALRA <- tsne_out$Y

umap_out <- umap(tsne_matrix)
umap_ALRA <- umap_out$layout
Sys.time()

#################################################### [8] ####################################################
method = "scGNGI"
method_scGNGI = "scGNGI"
file_name = "scGNGI_iter1_r3.csv"

# missing_impute_data_with_methods [2.0]
df_method <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,"_",method, "/scRNA_data3_group1_3_proprecessing_missing",i,"_",file_name))

rownames(df_method) <- df_method[,1]
df_method <- df_method[,-1]

# mask missing_zero_space 1
df_mask <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,"_mask.csv"))
rownames(df_mask) <- df_mask[,1]
df_mask <- df_mask[,-1]
df_mask <- as.matrix(df_mask)

df <- df_raw[,] + df_method*df_mask

df <- t(df)
df <- as.data.frame(df)
class <- substring(rownames(df), 1,6) # scRNA_data1_group1_5 // scRNA_data2_group1_4
label <- vector()
for (j in c(1:length(class)) ) {
  label[j] <- class[[j]][1]
}
df$label <- label

# 加载tsne包
library(Rtsne)
library(umap)
colors = rainbow(length(unique(df$label)))
names(colors) = unique(df$label)

# 使用tsne或者Rtsne函数进行tSNE降维分析
set.seed(42)
tsne_unique <- unique(df) # Remove duplicates
tsne_matrix <- as.matrix(tsne_unique[,1:(length(df)-1)])

Sys.time()
tsne_out <- Rtsne(tsne_matrix,dims=2,perplexity=30,theta=0.0,pca = FALSE) # Run TSNE
tsne_scGNGI <- tsne_out$Y

umap_out <- umap(tsne_matrix)
umap_scGNGI <- umap_out$layout
Sys.time()

save(tsne_raw, tsne_raw_missing, tsne_scImpute, tsne_SAVER, tsne_MAGIC, tsne_ALRA, tsne_DrImpute, tsne_scGNGI, umap_raw, umap_raw_missing, umap_scImpute, umap_SAVER, umap_MAGIC, umap_ALRA, umap_DrImpute, umap_scGNGI, file = paste0("./results/visualization/scRNA_data3_group1_3/tsne_umap_results_35.summary"))

#################################################### [7] ####################################################
setwd("~/r_workplace/8.genes_imputation/data")
load(paste0("./results/visualization/scRNA_data3_group1_3/tsne_umap_results_35.summary"))

#### plots tsne####
# 分割屏幕
split.screen(c(2, 4))

# 绘制第一个图
screen(1)
plot(tsne_raw_missing,col=colors[df$label],pch=16,
     xlab = "t-SNE1",ylab = "t-SNE2",main = method_raw_missing,cex.axis=0.93,cex.lab=0.93,cex.main=0.93)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 绘制第二个图
screen(2)
plot(tsne_scImpute,col=colors[df$label],pch=16,
     xlab = "t-SNE1",ylab = "t-SNE2",main = method_scImpute,cex.axis=0.93,cex.lab=0.93,cex.main=0.93)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 绘制第三个图
screen(3)
plot(tsne_SAVER,col=colors[df$label],pch=16,
     xlab = "t-SNE1",ylab = "t-SNE2",main = method_SAVER,cex.axis=0.93,cex.lab=0.93,cex.main=0.93)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 绘制第四个图
screen(4)
plot(tsne_MAGIC,col=colors[df$label],pch=16,
     xlab = "t-SNE1",ylab = "t-SNE2",main = method_MAGIC,cex.axis=0.93,cex.lab=0.93,cex.main=0.93)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 绘制第五个图
screen(5)
plot(tsne_ALRA,col=colors[df$label],pch=16,
     xlab = "t-SNE1",ylab = "t-SNE2",main = method_ALRA,cex.axis=0.93,cex.lab=0.93,cex.main=0.93)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 绘制第六个图
screen(6)
plot(tsne_DrImpute,col=colors[df$label],pch=16,
     xlab = "t-SNE1",ylab = "t-SNE2",main = method_DrImpute,cex.axis=0.93,cex.lab=0.93,cex.main=0.93)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 绘制第七个图
screen(7)
plot(tsne_scGNGI,col=colors[df$label],pch=16,
     xlab = "t-SNE1",ylab = "t-SNE2",main = method_scGNGI,cex.axis=0.93,cex.lab=0.93,cex.main=0.93)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 绘制第八个图
screen(8)
plot(tsne_raw,col=colors[df$label],pch=16,
     xlab = "t-SNE1",ylab = "t-SNE2",main = method_raw,cex.axis=0.93,cex.lab=0.93,cex.main=0.93)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 结束绘图
dev.off()

#### plots umap####
# 分割屏幕
split.screen(c(2, 4))

# 绘制第一个图
screen(1)
plot(umap_raw_missing,col=colors[df$label],pch=16,
     xlab = "UMAP_1",ylab = "UMAP_2",main = method_raw_missing,cex.axis=0.93,cex.lab=0.93,cex.main=0.93)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 绘制第二个图
screen(2)
plot(umap_scImpute,col=colors[df$label],pch=16,
     xlab = "UMAP_1",ylab = "UMAP_2",main = method_scImpute,cex.axis=0.93,cex.lab=0.93,cex.main=0.93)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 绘制第三个图
screen(3)
plot(umap_SAVER,col=colors[df$label],pch=16,
     xlab = "UMAP_1",ylab = "UMAP_2",main = method_SAVER,cex.axis=0.93,cex.lab=0.93,cex.main=0.93)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 绘制第四个图
screen(4)
plot(umap_MAGIC,col=colors[df$label],pch=16,
     xlab = "UMAP_1",ylab = "UMAP_2",main = method_MAGIC,cex.axis=0.93,cex.lab=0.93,cex.main=0.93)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 绘制第五个图
screen(5)
plot(umap_ALRA,col=colors[df$label],pch=16,
     xlab = "UMAP_1",ylab = "UMAP_2",main = method_ALRA,cex.axis=0.93,cex.lab=0.93,cex.main=0.93)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 绘制第六个图
screen(6)
plot(umap_DrImpute,col=colors[df$label],pch=16,
     xlab = "UMAP_1",ylab = "UMAP_2",main = method_DrImpute,cex.axis=0.93,cex.lab=0.93,cex.main=0.93)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 绘制第七个图
screen(7)
plot(umap_scGNGI,col=colors[df$label],pch=16,
     xlab = "UMAP_1",ylab = "UMAP_2",main = method_scGNGI,cex.axis=0.93,cex.lab=0.93,cex.main=0.93)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 绘制第八个图
screen(8)
plot(umap_raw,col=colors[df$label],pch=16,
     xlab = "UMAP_1",ylab = "UMAP_2",main = method_raw,cex.axis=0.93,cex.lab=0.93,cex.main=0.93)
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")

# 结束绘图
dev.off()

# 保存大小 width X height: pdf- 14 7,  eps- 1400 700








