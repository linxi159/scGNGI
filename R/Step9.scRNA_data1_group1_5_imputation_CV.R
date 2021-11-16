
###############################################
##### Coefficient of Variation(CV) 的分析 #####
###############################################

#################################################### [0] ####################################################
# load data
rm(list=ls())
setwd("~/r_workplace/8.genes_imputation/data")
# raw_data [1]
df_raw <- read.csv("./scRNA_data1_group1_5/scRNA_data1_group1_5_proprecessing.csv") 
rownames(df_raw) <- df_raw[,1]
df_raw <- df_raw[,-1]

##[1.1][2.1]
method1 = "scImpute"
file_name1 = "scImputescimpute_count.csv"

method2 = "SAVER"
file_name2 = "SAVER_estimate.csv"

method3 = "DrImpute"
file_name3 = "DrImpute.csv"

method4 = "ALRA"
file_name4 = "ALRA.csv"

method5 = "DeepImpute"
file_name5 = "DeepImpute.csv"

method6 = "scGNGI"
file_name6 = "scGNGI_iter1_r5.csv"


# raw_impute_data_with_methods [1.1]
df_method1 <- read.csv(paste0("./scRNA_data1_group1_5/scRNA_data1_group1_5_proprecessing_", method1, "/scRNA_data1_group1_5_proprecessing_",file_name1))
rownames(df_method1) <- df_method1[,1]
df_method1 <- df_method1[,-1]

df_method2 <- read.csv(paste0("./scRNA_data1_group1_5/scRNA_data1_group1_5_proprecessing_", method2, "/scRNA_data1_group1_5_proprecessing_",file_name2))
rownames(df_method2) <- df_method2[,1]
df_method2 <- df_method2[,-1]

df_method3 <- read.csv(paste0("./scRNA_data1_group1_5/scRNA_data1_group1_5_proprecessing_", method3, "/scRNA_data1_group1_5_proprecessing_",file_name3))
rownames(df_method3) <- df_method3[,1]
df_method3 <- df_method3[,-1]

df_method4 <- read.csv(paste0("./scRNA_data1_group1_5/scRNA_data1_group1_5_proprecessing_", method4, "/scRNA_data1_group1_5_proprecessing_",file_name4))
rownames(df_method4) <- df_method4[,1]
df_method4 <- df_method4[,-1]

df_method5 <- read.csv(paste0("./scRNA_data1_group1_5/scRNA_data1_group1_5_proprecessing_", method5, "/scRNA_data1_group1_5_proprecessing_",file_name5))
rownames(df_method5) <- df_method5[,1]
df_method5 <- df_method5[,-1]

df_method6 <- read.csv(paste0("./scRNA_data1_group1_5/scRNA_data1_group1_5_proprecessing_", method6, "/scRNA_data1_group1_5_proprecessing_",file_name6))
rownames(df_method6) <- df_method6[,1]
df_method6 <- df_method6[,-1]

log.scimpute.res<- log2(df_method1+1)
log.saver.res <- log2(df_method2+1)
log.drimpute.res <- log2(df_method4+1)
log.alra.res <- log2(df_method5+1)
log.deepimpute.res <- log2(df_method6+1)
log.scGNGI.res <- log2(df_method7+1)

#save(log.scimpute.res, log.saver.res, log.magic.res, log.drimpute.res,  log.alra.res, log.deepimpute.res, log.scGNGI.res, file = "./scRNA_data1_group1_5/scRNA_data1_group1_5_proprecessing_impute_res_summary")
save(log.scimpute.res, log.saver.res, log.drimpute.res,  log.alra.res, log.deepimpute.res, log.scGNGI.res, file = "./scRNA_data1_group1_5/scRNA_data1_group1_5_proprecessing_impute_res_summary")

#################################################### [1] ####################################################

load("./scRNA_data1_group1_5/scRNA_data1_group1_5_proprecessing_impute_res_summary")

CV.all <- function(data){
  CV <- function(mean, sd){
    abs(sd/mean)
  }
  CV.per.gene <- apply(data, 1, function(x){
    CV(mean(x), sd(x)) 
  })
  CV.per.gene
}	

scimpute.CV <- CV.all(log.scimpute.res)
saver.CV <- CV.all(log.saver.res)
drimpute.CV <- CV.all(log.drimpute.res)
alra.CV <- CV.all(log.alra.res)
deepimpute.CV <- CV.all(log.deepimpute.res)
scGNGI.CV <- CV.all(log.scGNGI.res)

celltype <- substring(colnames(df_raw), 1,6) # scRNA_data1_group1_5 // scRNA_data2_group1_4
label <- vector()
for (j in c(1:length(celltype)) ) {
  label[j] <- celltype[[j]][1]
}
celltype <- label

# raw_data df_raw
logxx <- apply(df_raw, 2, function(y){log2(y + 1)})
logxx[df_raw==0] <- 0

CV.nonzero <- function(data){
  CV <- function(mean, sd){
    abs(sd/mean)
  }
  CV.per.gene <- apply(data, 1, function(x){
    CV(mean(x[x!=0]), sd(x[x!=0])) 
  })
  CV.per.gene
}	
raw.CV <- CV.nonzero(logxx)

cell.category <- unique(celltype)
set.seed(5)
show.3000 <- sample(1:length(raw.CV), 3000)

library(ggplot2)
library(easyGgplot2)

for(i in 1:5){
  flag <- which(celltype%in%cell.category[i])
  raw.data <- logxx[, flag]
  zero.num <- apply(raw.data, 1, function(x){
    length(x[x==0])
  })
  zero.rate <- round(zero.num/ncol(raw.data), 2)*100
  dropout.rate <- apply(raw.data, 1, function(x){
    round(mean(x[x!=0]), 2)
  })
  
  scimpute.CV <- CV.all(log.scimpute.res[, flag])[show.3000]
  saver.CV <- CV.all(log.saver.res[, flag])[show.3000]
  drimpute.CV <- CV.all(log.drimpute.res[, flag])[show.3000]
  alra.CV <- CV.all(log.alra.res[, flag])[show.3000]
  deepimpute.CV <- CV.all(log.deepimpute.res[, flag])[show.3000]
  scGNGI.CV <- CV.all(log.scGNGI.res[, flag])[show.3000]
  raw.CV <- CV.nonzero(logxx[, flag])[show.3000]
  zero.rate.selected <- zero.rate[show.3000]
  dropout.rate.selected <- dropout.rate[show.3000]
  
  df <- data.frame(value = c(scimpute.CV, saver.CV, drimpute.CV, alra.CV, deepimpute.CV, scGNGI.CV), method=c(rep("scImpute", 3000), rep("SAVER", 3000), rep("DrImpute", 3000),  rep("ALRA", 3000), rep("DeepImpute", 3000),rep("scGNGI", 3000)), compr = rep("Without-imputation", 3000*6), raw = rep(raw.CV, 6), zero = rep(zero.rate.selected, 6), dropout = rep(dropout.rate.selected, 6))

  gg <- ggplot2.scatterplot(data = df, xName = 'raw', yName = 'value', size = 0.3, backgroundColor = "white", xtitle="CV (Before Imputation)", 
                            ytitle="CV (After Imputation)", mainTitle = cell.category[i], removePanelGrid=TRUE,   removePanelBorder=FALSE, showLegend=TRUE, 
                            legendTitle = "Percentage \n of Zero", legendTitleFont = c(10, "bold", "black"), legendTextFont = c(10, "bold", "black"), 
                            mainTitleFont = c(10, "bold", "black"), xtitleFont = c(10, "bold", "black"),  ytitleFont = c(10, "bold", "black"), 
                            xTickLabelFont = c(10, "bold", "white"), yTickLabelFont = c(10, "bold", "black"), facetingFont = c(10, "bold", "black"), 
                            facetingRect =  list(background = NULL, lineType = NULL, lineColor = NULL, lineSize = NULL)) + geom_point(aes(colour = zero), size = 0.3) + ylim(0, 5)  + xlim(0, 5) + facet_wrap(~method, ncol=3)  + theme(strip.text.x = element_text(size = 10, colour = "black", face = "bold")) + geom_abline(col = "red", linetype = "dashed", size=1) 
  ggsave(gg, file = paste0("./results/CV/scRNA_data1_group1_5_proprecessing_",  cell.category[i], "_CV_zero_comparison.eps"), width = 12, height = 7)
  

}




