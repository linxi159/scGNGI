
####@@@@ missing 10 5 2 @@@@####
rm(list=ls())
#### summary ####
i = 10 # 10 5 2
setwd("~/r_workplace/8.genes_imputation/data")
# 1. scRNA_data1_group1_5_proprecessing.csv 原始的预处理数据
df0 <- read.csv("./scRNA_data1_group1_5/scRNA_data1_group1_5_proprecessing.csv")
rownames(df0) <- df0[,1]
df0 <- df0[,-1]

# 2. scRNA_data1_group1_5_proprecessing_missingi_mask.csv 掩码位置
df0_missing10.mask <- read.csv(paste0("./scRNA_data1_group1_5/scRNA_data1_group1_5_proprecessing_missing",i,"_mask.csv"))
rownames(df0_missing10.mask) <- df0_missing10.mask[,1]
df0_missing10.mask <- df0_missing10.mask[,-1]

# 3. scRNA_data1_group1_5_proprecessing_missingi_methods/...csv 填充后数据
df0_missing10.scimpute_imputed <- read.csv(paste0("./scRNA_data1_group1_5/scRNA_data1_group1_5_proprecessing_missing",i,"_scImpute/scRNA_data1_group1_5_proprecessing_missing",i,"_scImputescimpute_count.csv"))
df0_missing10.saver_imputed <- read.csv(paste0("./scRNA_data1_group1_5/scRNA_data1_group1_5_proprecessing_missing",i,"_SAVER/scRNA_data1_group1_5_proprecessing_missing",i,"_SAVER_estimate.csv"))
df0_missing10.magic_imputed <- read.csv(paste0("./scRNA_data1_group1_5/scRNA_data1_group1_5_proprecessing_missing",i,"_MAGIC/scRNA_data1_group1_5_proprecessing_missing",i,"_MAGIC.csv"))
df0_missing10.drimpute_imputed <- read.csv(paste0("./scRNA_data1_group1_5/scRNA_data1_group1_5_proprecessing_missing",i,"_DrImpute/scRNA_data1_group1_5_proprecessing_missing",i,"_DrImpute.csv"))
df0_missing10.alra_imputed <- read.csv(paste0("./scRNA_data1_group1_5/scRNA_data1_group1_5_proprecessing_missing",i,"_ALRA/scRNA_data1_group1_5_proprecessing_missing",i,"_ALRA.csv"))
df0_missing10.deepimpute_imputed <- read.csv(paste0("./scRNA_data1_group1_5/scRNA_data1_group1_5_proprecessing_missing",i,"_DeepImpute/scRNA_data1_group1_5_proprecessing_missing",i,"_DeepImpute.csv"))
df0_missing10.scGNGI_imputed <- read.csv(paste0("./scRNA_data1_group1_5/scRNA_data1_group1_5_proprecessing_missing",i,"_scGNGI/scRNA_data1_group1_5_proprecessing_missing",i,"_scGNGI_iter1_r5.csv"))

rownames(df0_missing10.scimpute_imputed) <- df0_missing10.scimpute_imputed[,1]
df0_missing10.scimpute_imputed <- df0_missing10.scimpute_imputed[,-1]
rownames(df0_missing10.saver_imputed) <- df0_missing10.saver_imputed[,1]
df0_missing10.saver_imputed <- df0_missing10.saver_imputed[,-1]
rownames(df0_missing10.magic_imputed) <- df0_missing10.magic_imputed[,1]
df0_missing10.magic_imputed <- df0_missing10.magic_imputed[,-1]
rownames(df0_missing10.drimpute_imputed) <- df0_missing10.drimpute_imputed[,1]
df0_missing10.drimpute_imputed <- df0_missing10.drimpute_imputed[,-1]
rownames(df0_missing10.alra_imputed) <- df0_missing10.alra_imputed[,1]
df0_missing10.alra_imputed <- df0_missing10.alra_imputed[,-1]
rownames(df0_missing10.deepimpute_imputed) <- df0_missing10.deepimpute_imputed[,1]
df0_missing10.deepimpute_imputed <- df0_missing10.deepimpute_imputed[,-1]
rownames(df0_missing10.scGNGI_imputed) <- df0_missing10.scGNGI_imputed[,1]
df0_missing10.scGNGI_imputed <- df0_missing10.scGNGI_imputed[,-1]

# 提取插补的数值
# a=as.matrix(df0[1:17,1:4])
# b=as.matrix(df0_missing10.mask[1:17,1:4])
# c=a*b
# d=as.numeric(c) 
# to=d[d!=0] 按列提取分别存储

# 原始数据中对应missing的部分 
topredict <- as.numeric(as.matrix(df0)) 
missing <-  as.numeric(as.matrix(df0_missing10.mask))
topredict <- topredict[missing!=0]
# 插补数据中对应missing的部分 
scimpute.pred <- as.numeric(as.matrix(df0_missing10.scimpute_imputed))
scimpute.pred <- scimpute.pred[missing!=0]
saver.pred <- as.numeric(as.matrix(df0_missing10.saver_imputed))
saver.pred <- saver.pred[missing!=0]
magic.pred <- as.numeric(as.matrix(df0_missing10.magic_imputed))
magic.pred <- magic.pred[missing!=0]
drimpute.pred <- as.numeric(as.matrix(df0_missing10.drimpute_imputed))
drimpute.pred <- drimpute.pred[missing!=0]
alra.pred <- as.numeric(as.matrix(df0_missing10.alra_imputed))
alra.pred <- alra.pred[missing!=0]
deepimpute.pred <- as.numeric(as.matrix(df0_missing10.deepimpute_imputed))
deepimpute.pred <- deepimpute.pred[missing!=0]
scGNGI.pred <- as.numeric(as.matrix(df0_missing10.scGNGI_imputed))
scGNGI.pred <- scGNGI.pred[missing!=0]

# save
save(topredict, scimpute.pred, saver.pred, magic.pred, drimpute.pred, alra.pred,deepimpute.pred, scGNGI.pred, file = paste0("./scRNA_data1_group1_5/results_missing",i,"/scRNA_data1_group1_5_proprecessing_missing",i,"_imputation_log2.summary"))

#### evaluate ####
setwd("~/r_workplace/8.genes_imputation/data")
# absolute.loss
all.res <- NULL
load(paste0("./scRNA_data1_group1_5/results_missing",i,"/scRNA_data1_group1_5_proprecessing_missing",i,"_imputation_log2.summary"))

#数据规范化 [log2(x+1)] log10 log
# log.topredict <- log(topredict+1) // log.topredict <- log10(topredict+1)

log.topredict <- log2(topredict+1)
log.scimpute.pred <- log2(scimpute.pred+1)
log.saver.pred <- log2(saver.pred+1)
log.magic.pred <- log2(magic.pred+1)
log.drimpute.pred <- log2(drimpute.pred+1)
log.alra.pred <- log2(alra.pred+1)
log.deepimpute.pred <- log2(deepimpute.pred+1)
log.scGNGI.pred <- log2(scGNGI.pred+1)

aa <- c(mean(abs(log.topredict - log.scimpute.pred)),
        mean(abs(log.topredict - log.saver.pred)),
        mean(abs(log.topredict - log.magic.pred)),
        mean(abs(log.topredict - log.drimpute.pred)),
        mean(abs(log.topredict - log.alra.pred)),
        mean(abs(log.topredict - log.deepimpute.pred)),
        mean(abs(log.topredict - log.scGNGI.pred)))

all.res <- rbind(all.res, aa)

absolute.loss.mean.summary <- data.frame(dataset = rep("Simutation Data1", nrow(all.res)), 
                                         scenario = rep(paste0("missing",i), nrow(all.res)),
                                         scimpute = all.res[, 1],
                                         saver = all.res[, 2],
                                         magic = all.res[, 3],
                                         drimpute = all.res[, 4],
                                         alra  = all.res[, 5],
                                         deepimpute = all.res[, 6],
                                         scGNGI = all.res[, 7]) 

# square.loss 的均值
all.res <- NULL
load(paste0("./scRNA_data1_group1_5/results_missing",i,"/scRNA_data1_group1_5_proprecessing_missing",i,"_imputation_log2.summary"))

#数据规范化 [log2(x+1)] log10 log
# log.topredict <- log(topredict+1) // log.topredict <- log10(topredict+1)

log.topredict <- log2(topredict+1)
log.scimpute.pred <- log2(scimpute.pred+1)
log.saver.pred <- log2(saver.pred+1)
log.magic.pred <- log2(magic.pred+1)
log.drimpute.pred <- log2(drimpute.pred+1)
log.alra.pred <- log2(alra.pred+1)
log.deepimpute.pred <- log2(deepimpute.pred+1)
log.scGNGI.pred <- log2(scGNGI.pred+1)

aa <- c(mean((log.topredict - log.scimpute.pred)^2),
        mean((log.topredict - log.saver.pred)^2),
        mean((log.topredict - log.magic.pred)^2),
        mean((log.topredict - log.drimpute.pred)^2),
        mean((log.topredict - log.alra.pred)^2),
        mean((log.topredict - log.deepimpute.pred)^2),
        mean((log.topredict - log.scGNGI.pred)^2))

all.res <- rbind(all.res, aa)

square.loss.mean.summary <- data.frame(dataset = rep("Simutation Data1", nrow(all.res)), 
                                       scenario = rep(paste0("missing,",i), nrow(all.res)),
                                       scimpute = all.res[, 1],
                                       saver = all.res[, 2],
                                       magic = all.res[, 3],
                                       drimpute = all.res[, 4],
                                       alra  = all.res[, 5],
                                       deepimpute = all.res[, 6],
                                       scGNGI = all.res[, 7]) 
# correlation
all.res <- NULL
load(paste0("./scRNA_data1_group1_5/results_missing",i,"/scRNA_data1_group1_5_proprecessing_missing",i,"_imputation_log2.summary"))

#数据规范化 [log2(x+1)] log10 log
# log.topredict <- log(topredict+1) // log.topredict <- log10(topredict+1)

log.topredict <- log2(topredict+1)
log.scimpute.pred <- log2(scimpute.pred+1)
log.saver.pred <- log2(saver.pred+1)
log.magic.pred <- log2(magic.pred+1)
log.drimpute.pred <- log2(drimpute.pred+1)
log.alra.pred <- log2(alra.pred+1)
log.deepimpute.pred <- log2(deepimpute.pred+1)
log.scGNGI.pred <- log2(scGNGI.pred+1)

aa <- c(cor(log.topredict, log.scimpute.pred),
        cor(log.topredict, log.saver.pred),
        cor(log.topredict, log.magic.pred),
        cor(log.topredict, log.drimpute.pred),
        cor(log.topredict, log.alra.pred),
        cor(log.topredict, log.deepimpute.pred),
        cor(log.topredict, log.scGNGI.pred))

all.res <- rbind(all.res, aa)

correlation.mean.summary <- data.frame(dataset = rep("Simutation Data1", nrow(all.res)), 
                                       scenario = rep(paste0("missing,",i), nrow(all.res)),
                                       scimpute = all.res[, 1],
                                       saver = all.res[, 2],
                                       magic = all.res[, 3],
                                       drimpute = all.res[, 4],
                                       alra  = all.res[, 5],
                                       deepimpute = all.res[, 6],
                                       scGNGI = all.res[, 7]) 

# F范数1-对应元素的平方和再开方, F范数(缺失部分的插补前后的差值)/F范数(原始数据)
all.res <- NULL
load(paste0("./scRNA_data1_group1_5/results_missing",i,"/scRNA_data1_group1_5_proprecessing_missing",i,"_imputation_log2.summary"))

#数据规范化 [log2(x+1)] log10 log
# log.topredict <- log(topredict+1) // log.topredict <- log10(topredict+1)

log.topredict <- log2(topredict+1)
log.scimpute.pred <- log2(scimpute.pred+1)
log.saver.pred <- log2(saver.pred+1)
log.magic.pred <- log2(magic.pred+1)
log.drimpute.pred <- log2(drimpute.pred+1)
log.alra.pred <- log2(alra.pred+1)
log.deepimpute.pred <- log2(deepimpute.pred+1)
log.scGNGI.pred <- log2(scGNGI.pred+1)

aa <- c(cor(log.topredict, log.scimpute.pred),
        cor(log.topredict, log.saver.pred),
        cor(log.topredict, log.magic.pred),
        cor(log.topredict, log.drimpute.pred),
        cor(log.topredict, log.alra.pred),
        cor(log.topredict, log.deepimpute.pred),
        cor(log.topredict, log.scGNGI.pred))

df0 <- read.csv("./scRNA_data1_group1_5/scRNA_data1_group1_5_proprecessing.csv")
rownames(df0) <- df0[,1]
df0 <- df0[,-1]
df0_ <- as.numeric(as.matrix(df0)) 
log.df0 <- log2(df0_+1)

aa <- c(sqrt(sum((log.topredict - log.scimpute.pred)^2))/sqrt(sum((log.df0)^2)),
        sqrt(sum((log.topredict - log.saver.pred)^2))/sqrt(sum((log.df0)^2)),
        sqrt(sum((log.topredict - log.magic.pred)^2))/sqrt(sum((log.df0)^2)),
        sqrt(sum((log.topredict - log.drimpute.pred)^2))/sqrt(sum((log.df0)^2)),
        sqrt(sum((log.topredict - log.alra.pred)^2))/sqrt(sum((log.df0)^2)),
        sqrt(sum((log.topredict - log.deepimpute.pred)^2))/sqrt(sum((log.df0)^2)),
        sqrt(sum((log.topredict - log.scGNGI.pred)^2))/sqrt(sum((log.df0)^2)))

all.res <- rbind(all.res, aa)

f.loss.mean.summary <- data.frame(dataset = rep("Simutation Data1", nrow(all.res)), 
                                  scenario = rep(paste0("missing,",i), nrow(all.res)),
                                  scimpute = all.res[, 1],
                                  saver = all.res[, 2],
                                  magic = all.res[, 3],
                                  drimpute = all.res[, 4],
                                  alra  = all.res[, 5],
                                  deepimpute = all.res[, 6],
                                  scGNGI = all.res[, 7]) 

# save
save(f.loss.mean.summary, correlation.mean.summary, square.loss.mean.summary, absolute.loss.mean.summary, file = paste0("./scRNA_data1_group1_5/results_missing",i,"/scRNA_data1_group1_5_proprecessing_missing",i,"_evaluate_log2.summary"))


#### plots ####


