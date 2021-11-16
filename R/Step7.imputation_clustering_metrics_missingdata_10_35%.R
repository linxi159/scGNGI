
###############################################
############## Clustering  Examples ###########
###############################################

#################################################### [0] ####################################################
# load data
rm(list=ls())
setwd("~/r_workplace/8.genes_imputation/data")

i = 10 # 10 35
# raw_data [1]
# df_raw <- read.csv("./GSE75748/GSE75748_sc_cell_type_ec_proprecessing.csv")
# raw_data_missing2510 [2]
df_raw <- read.csv(paste0("./GSE75748/GSE75748_sc_cell_type_ec_proprecessing_missing",i,".csv"))

# i = 10 # 10 35
# # raw_data [1]
# #df_raw <- read.csv("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing.csv") 
# # raw_data_missing2510 [2]
# df_raw <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,".csv")) 

# i = 10 # 10 
# # raw_data [1]
# #df_raw <- read.csv("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing.csv") 
# # raw_data_missing2510 [2]
# df_raw <- read.csv(paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missing",i,".csv")) 

# i = 10 # 10 
# raw_data [1]
#df_raw <- read.csv("./scRNA_data1_group1_5/scRNA_data1_group1_5_proprecessing.csv") 
# raw_data_missing2510 [2]
# df_raw <- read.csv(paste0("./scRNA_data1_group1_5/scRNA_data1_group1_5_proprecessing_missing",i,".csv")) 

rownames(df_raw) <- df_raw[,1]
df_raw <- df_raw[,-1]

##[1.1][2.1]
# method = "scImpute"
# file_name = "scImputescimpute_count.csv"
# 
# method = "SAVER"
# file_name = "SAVER_estimate.csv"
# 
# method = "MAGIC"
# file_name = "MAGIC.csv"
# 
# method = "DrImpute"
# file_name = "DrImpute.csv"
# 
# method = "ALRA"
# file_name = "ALRA.csv"
# 
method = "scGNGI"
file_name = "scGNGI_r7.csv"

# raw_impute_data_with_methods [1.1]
#df_method <- read.csv(paste0("./GSE75748/GSE75748_sc_cell_type_ec_proprecessing_", method, "/GSE75748_sc_cell_type_ec_proprecessing_",file_name))
# missing_impute_data_with_methods [2.1]
df_method <- read.csv(paste0("./GSE75748/GSE75748_sc_cell_type_ec_proprecessing_missing",i,"_",method, "/GSE75748_sc_cell_type_ec_proprecessing_missing",i,"_",file_name))

# missing_impute_data_with_methods [2.1]
# df_method <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,"_",method, "/scRNA_data3_group1_3_proprecessing_missing",i,"_",file_name))

# missing_impute_data_with_methods [2.1]
# df_method <- read.csv(paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missing",i,"_",method, "/scRNA_data2_group1_4_proprecessing_missing",i,"_",file_name))

# missing_impute_data_with_methods [2.1]
# df_method <- read.csv(paste0("./scRNA_data1_group1_5/scRNA_data1_group1_5_proprecessing_missing",i,"_",method, "/scRNA_data1_group1_5_proprecessing_missing",i,"_",file_name))

rownames(df_method) <- df_method[,1]
df_method <- df_method[,-1]

# mask missing_zero_space 1
df_mask <- read.csv(paste0("./GSE75748/GSE75748_sc_cell_type_ec_proprecessing_missing",i,"_mask.csv"))
# df_mask <- read.csv(paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing",i,"_mask.csv"))
# df_mask <- read.csv(paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missing",i,"_mask.csv"))
# df_mask <- read.csv(paste0("./scRNA_data1_group1_5/scRNA_data1_group1_5_proprecessing_missing",i,"_mask.csv"))
rownames(df_mask) <- df_mask[,1]
df_mask <- df_mask[,-1]
df_mask <- as.matrix(df_mask)

df <- df_raw[,] + df_method*df_mask
#df <- df_method

df <- t(df)
df <- as.data.frame(df)
class <- strsplit(rownames(df), "_") # GSE75748
#class <- substring(rownames(df), 1,6) # scRNA_data3_group1_3 
label <- vector()
for (j in c(1:length(class)) ) {
  label[j] <- class[[j]][1]
}


#################################################### [1] ####################################################
# NMI 计算
library(infotheo)
library(Rtsne)
library(umap)
set.seed(42)

# 使用Rtsne或者Umap函数进行降维
set.seed(42)
tsne_unique <- unique(df) # Remove duplicates
tsne_matrix <- as.matrix(tsne_unique[,1:length(df)])

Sys.time()
tsne_out <- Rtsne(tsne_matrix,dims=2,perplexity=30,theta=0.0,pca = FALSE) # Run TSNE
tsne_df <- tsne_out$Y

umap_out <- umap(tsne_matrix)
umap_df <- umap_out$layout
Sys.time()

kc <- kmeans(df, 7)
kc_tsne <- kmeans(tsne_df, 7);  #分类模型训练
kc_umap <- kmeans(umap_df, 7);  #分类模型训练

label_prediction=kc$cluster
label_prediction_tsne=kc_tsne$cluster
label_prediction_umap=kc_umap$cluster
for (i in 1:length(label)) {
  if (label[i] == "H1"){
    label[i] =1
  }
  if (label[i] == "H9"){
    label[i] =2
  }
  if (label[i] == "DEC"){
    label[i] =3
  }
  if (label[i] == "EC"){
    label[i] =4
  }
  if (label[i] == "HFF"){
    label[i] =5
  }
  if (label[i] == "NPC"){
    label[i] =6
  }
  if (label[i] == "TB"){
    label[i] =7
  }
}
# for (i in 1:length(label)) {
#   if (label[i] == "Group1"){
#     label[i] =1
#   }
#   if (label[i] == "Group2"){
#     label[i] =2
#   }
#   if (label[i] == "Group3"){
#     label[i] =3
#   }
# }
# for (i in 1:length(label)) {
#   if (label[i] == "Group1"){
#     label[i] =1
#   }
#   if (label[i] == "Group2"){
#     label[i] =2
#   }
#   if (label[i] == "Group3"){
#     label[i] =3
#   }
#   if (label[i] == "Group4"){
#     label[i] =4
#   }
# }
# for (i in 1:length(label)) {
#   if (label[i] == "Group1"){
#     label[i] =1
#   }
#   if (label[i] == "Group2"){
#     label[i] =2
#   }
#   if (label[i] == "Group3"){
#     label[i] =3
#   }
#   if (label[i] == "Group4"){
#     label[i] =4
#   }
#   if (label[i] == "Group5"){
#     label[i] =5
#   }
# }
label = as.numeric(label)
label_prediction = as.numeric(label_prediction)
label_prediction_tsne = as.numeric(label_prediction_tsne)
label_prediction_umap = as.numeric(label_prediction_umap)

nmi = 2 * mutinformation(label,label_prediction) /(entropy(label)+entropy(label_prediction))
nmi_tsne = 2 * mutinformation(label,label_prediction_tsne) /(entropy(label)+entropy(label_prediction_tsne))
nmi_umap = 2 * mutinformation(label,label_prediction_umap) /(entropy(label)+entropy(label_prediction_umap))
nmi
nmi_tsne
nmi_umap 

#################################################### [2] ####################################################
# plot figure 画图展示
# scRNA_data3_missing10
library(ggplot2)
library(gridExtra)
Method <- c("scImpute", "SAVER", "MAGIC", "DrImpute","ALRA","scGNGI")
nmi <- c(0.489, 0.476, 0.342, 0.476, 0.500,0.501)
scRNA_data3_missing10 <- data.frame(Method, nmi)

pl1 <- ggplot(data=scRNA_data3_missing10,aes(x=Method,y=nmi,fill=Method),fill=Method) +
  geom_bar(stat="identity") + theme(panel.background=element_rect(fill='transparent',color ="black"),axis.text.x = element_blank(),panel.grid =element_blank()) +
  labs(x = "Imputation Methods", y= "NMI") +
  geom_text(aes(label=nmi, y=nmi+0.01), position=position_dodge(0.9), vjust=0) +
  ggtitle("(c) Simulated Data 3 with Missing 10%") + theme(plot.title = element_text(size = 10.5, face = "bold"))
                                 
# scRNA_data3_missing35
Method <- c("scImpute", "SAVER", "MAGIC", "DrImpute","ALRA","scGNGI")
nmi <- c(0.224, 0.138, 0.256, 0.258, 0.254,0.285)
scRNA_data3_missing35 <- data.frame(Method, nmi)

pl2 <- ggplot(data=scRNA_data3_missing35,aes(x=Method,y=nmi,fill=Method),fill=Method) +
  geom_bar(stat="identity") + theme(panel.background=element_rect(fill='transparent',color ="black"),axis.text.x = element_blank(),panel.grid =element_blank()) +
  labs(x = "Imputation Methods", y= "NMI") +
  geom_text(aes(label=nmi, y=nmi+0.01), position=position_dodge(0.9), vjust=0) +
  ggtitle("(d) Simulated Data 3 with Missing 35%") + theme(plot.title = element_text(size = 10.5, face = "bold"))

# scRNA_data1_missing10
Method <- c("scImpute", "SAVER", "MAGIC", "DrImpute","ALRA","scGNGI")
nmi <- c( 0.99, 0.981, 0.987, 0.978, 0.851,0.987)
scRNA_data2_missing10 <- data.frame(Method, nmi)

pl3 <- ggplot(data=scRNA_data2_missing10,aes(x=Method,y=nmi,fill=Method),fill=Method) +
  geom_bar(stat="identity") + theme(panel.background=element_rect(fill='transparent',color ="black"),axis.text.x = element_blank(),panel.grid =element_blank()) +
  labs(x = "Imputation Methods", y= "NMI") +
  geom_text(aes(label=nmi, y=nmi+0.01), position=position_dodge(0.9), vjust=0) +
  ggtitle("(b) Simulated Data 1 with Missing 10%") + theme(plot.title = element_text(size = 10.5, face = "bold"))

# GSE75748_missing10
Method <- c("scImpute", "SAVER", "MAGIC", "DrImpute","ALRA","scGNGI")
nmi <- c(0.58, 0.57, 0.563, 0.556, 0.574, 0.576)
real_data_missing10 <- data.frame(Method, nmi)

pl4 <- ggplot(data=real_data_missing10,aes(x=Method,y=nmi,fill=Method),fill=Method) +
  geom_bar(stat="identity") + theme(panel.background=element_rect(fill='transparent',color ="black"),axis.text.x = element_blank(),panel.grid =element_blank()) +
  labs(x = "Imputation Methods", y= "NMI") +
  geom_text(aes(label=nmi, y=nmi+0.01), position=position_dodge(0.9), vjust=0) +
  ggtitle("(a) Cell Type with Missing 10%") + theme(plot.title = element_text(size = 10.5, face = "bold"))

grid.arrange(pl4,pl3,pl1,pl2,ncol=2)
# 保存大小 width X height: pdf- 12 9,  eps- 1200 900

