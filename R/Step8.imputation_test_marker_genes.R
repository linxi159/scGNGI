
##########################################################
##### Marker genes non-linear relationships analysis #####
##########################################################

#################################################### [0] ####################################################
# load data
rm(list=ls())
setwd("~/r_workplace/8.genes_imputation/data")

# raw_data [1]
df_raw <- read.csv("./GSE75748/GSE75748_sc_cell_type_ec_proprecessing.csv")
rownames(df_raw) <- df_raw[,1]
df_raw <- df_raw[,-1]

##[1.1][2.1]
method1 = "scImpute"
file_name1 = "scImputescimpute_count.csv"

method2 = "SAVER"
file_name2 = "SAVER_estimate.csv"

method3 = "MAGIC"
file_name3 = "MAGIC.csv"

method4 = "ALRA"
file_name4 = "ALRA.csv"

method5 = "DeepImpute"
file_name5 = "DeepImpute.csv"

method6 = "scGNGI"
file_name6 = "scGNGI_iter1_r7.csv"


# raw_impute_data_with_methods [1.1]
df_method1 <- read.csv(paste0("./GSE75748/GSE75748_sc_cell_type_ec_proprecessing_", method1, "/GSE75748_sc_cell_type_ec_proprecessing_",file_name1))
rownames(df_method1) <- df_method1[,1]
df_method1 <- df_method1[,-1]

df_method2 <- read.csv(paste0("./GSE75748/GSE75748_sc_cell_type_ec_proprecessing_", method2, "/GSE75748_sc_cell_type_ec_proprecessing_",file_name2))
rownames(df_method2) <- df_method2[,1]
df_method2 <- df_method2[,-1]

df_method3 <- read.csv(paste0("./GSE75748/GSE75748_sc_cell_type_ec_proprecessing_", method3, "/GSE75748_sc_cell_type_ec_proprecessing_",file_name3))
rownames(df_method3) <- df_method3[,1]
df_method3 <- df_method3[,-1]

df_method4 <- read.csv(paste0("./GSE75748/GSE75748_sc_cell_type_ec_proprecessing_", method4, "/GSE75748_sc_cell_type_ec_proprecessing_",file_name4))
rownames(df_method4) <- df_method4[,1]
df_method4 <- df_method4[,-1]

df_method5 <- read.csv(paste0("./GSE75748/GSE75748_sc_cell_type_ec_proprecessing_", method5, "/GSE75748_sc_cell_type_ec_proprecessing_",file_name5))
rownames(df_method5) <- df_method5[,1]
df_method5 <- df_method5[,-1]

df_method6 <- read.csv(paste0("./GSE75748/GSE75748_sc_cell_type_ec_proprecessing_", method6, "/GSE75748_sc_cell_type_ec_proprecessing_",file_name6))
rownames(df_method6) <- df_method6[,1]
df_method6 <- df_method6[,-1]


#################################################### [1] ####################################################
# 对于GSE75748数据的cell marker的gene pairs(基因对)的非线性相关性恢复的分析
cell_marker <- c("KLF4","NANOG", "SOX2","CD9","CDH11","EFNA2","PRSS50")

library(minerva)
# 检验非线性相关性(MIC度量)和P值(wilcox.test)
mine_raw<-c()
mine_raw_p<-c()
mine_method<-c()
mine_method_p<-c()
df_method = df_method6 # 6 2 5 

for (i in cell_marker) {
  for (j in cell_marker) {
    if( i != j ){
      x = scale(as.numeric(df_raw[i,]))
      y = scale(as.numeric(df_raw[j,]))
      a=mine(x,y)
      mine_raw = c(mine_raw, a$MIC)
      tmp=wilcox.test(x,y)
      mine_raw_p = c(mine_raw_p,tmp[["p.value"]] )
        
      x_ = scale(as.numeric(df_method[i,]))
      y_ = scale(as.numeric(df_method[j,]))
      b=mine(x_,y_)
      mine_method = c(mine_method, b$MIC )
      tmp_ = wilcox.test(x_,y_)
      mine_method_p = c(mine_method_p,tmp_[["p.value"]] )
      
    }
  }
}

# mine_raw
# mine_raw_p
# mine_method
# mine_method_p
# 
# mean(mine_raw)
# mean(mine_method)

tmp=wilcox.test(mine_raw,mine_method)
tmp[["p.value"]]

##################### [2] # plot #################
# NANOG - PRSS50
# EFNA2 - PRSS50
# SOX2 - PRSS50
# CDH11 - PRSS50

i = "CDH11" # NANOG EFNA2 SOX2 CDH11
j= "PRSS50"

# 归一化
normlize <- function(x){
  a = (x - min(x)) / (max(x) - min(x))
  return(a)
}

mic_p <- function(x_,y_){
  mine_raw<-c()
  mine_raw_p<-c()
  x = scale(as.numeric(x_))
  y = scale(as.numeric(y_))
  a=mine(x,y)
  mine_raw = c(mine_raw, a$MIC)
  tmp=wilcox.test(x,y)
  mine_raw_p = c(mine_raw_p,tmp[["p.value"]] )

  return (list(mine_raw, mine_raw_p))
}

mine_raw<-c()
mine_raw_p<-c()
mine_method<-c()
mine_method_p<-c()

# 分割屏幕
split.screen(c(2, 3))
# 绘制第一个图
screen(1)

library(stringr)
a = mic_p(df_raw[i,],df_raw[j,])
mine_raw=as.character(a[1])
mine_raw_p=as.character(a[2])
mine_raw = str_c(substr(mine_raw,start=1 ,stop=5),",")
mine_raw_p = str_c(substr(mine_raw_p,start=1 ,stop=5),"e-",strsplit(mine_raw_p,split = "-")[[1]][2])
cc=paste("MIC = ",mine_raw,"P.value = ",mine_raw_p)

xx = normlize(as.numeric(df_raw[i,]))
yy = normlize(as.numeric(df_raw[j,]))
name = "Raw"

plot(xx, yy, main=name, xlab=i, ylab=j, pch=19) #绘图
#legend("topleft", inset=0.001,c(cc),bty = "n",title.adj = 1)
legend(x=-0.83,y=1.24, c(cc),bty = "n")


# 绘制第二个图
screen(2)

df_method = df_method1 # 1 2 5 6 7
a = mic_p(df_method[i,],df_method[j,])
mine_method=as.character(a[1])
mine_method_p=as.character(a[2])
mine_method = str_c(substr(mine_method,start=1 ,stop=5),",")
mine_method_p = str_c(substr(mine_method_p,start=1 ,stop=5),"e-",strsplit(mine_method_p,split = "-")[[1]][2])
cc=paste("MIC = ",mine_method,"P.value = ",mine_method_p)

xx = normlize(as.numeric(df_method[i,]))
yy = normlize(as.numeric(df_method[j,]))
name = "scImpute"

plot(xx, yy, main=name, xlab=i, ylab=j, pch=19) #绘图
legend(x=-0.83,y=1.24, c(cc),bty = "n")

# 绘制第三个图
screen(3)

df_method = df_method2 # 1 2 5 6 7
a = mic_p(df_method[i,],df_method[j,])
mine_method=as.character(a[1])
mine_method_p=as.character(a[2])
mine_method = str_c(substr(mine_method,start=1 ,stop=5),",")
mine_method_p = str_c(substr(mine_method_p,start=1 ,stop=5),"e-",strsplit(mine_method_p,split = "-")[[1]][2])
cc=paste("MIC = ",mine_method,"P.value = ",mine_method_p)

xx = normlize(as.numeric(df_method[i,]))
yy = normlize(as.numeric(df_method[j,]))
name = "SAVER"

plot(xx, yy, main=name, xlab=i, ylab=j, pch=19) #绘图
legend(x=-0.83,y=1.24, c(cc),bty = "n")

# 绘制第四个图
screen(4)

df_method = df_method5 # 1 2 5 6 7
a = mic_p(df_method[i,],df_method[j,])
mine_method=as.character(a[1])
mine_method_p=as.character(a[2])
mine_method = str_c(substr(mine_method,start=1 ,stop=5),",")
mine_method_p = str_c(substr(mine_method_p,start=1 ,stop=5),"e-",strsplit(mine_method_p,split = "-")[[1]][2])
cc=paste("MIC = ",mine_method,"P.value = ",mine_method_p)

xx = normlize(as.numeric(df_method[i,]))
yy = normlize(as.numeric(df_method[j,]))
name = "ALRA"

plot(xx, yy, main=name, xlab=i, ylab=j, pch=19) #绘图
legend(x=-0.83,y=1.24, c(cc),bty = "n")

# 绘制第五个图
screen(5)

df_method = df_method6 # 1 2 5 6 7
a = mic_p(df_method[i,],df_method[j,])
mine_method=as.character(a[1])
mine_method_p=as.character(a[2])
mine_method = str_c(substr(mine_method,start=1 ,stop=5),",")
#mine_method_p = str_c(substr(mine_method_p,start=1 ,stop=5),"e-",strsplit(mine_method_p,split = "-")[[1]][2])
#mine_method_p = "8.916e-01"
#mine_method_p = "6.419e-01"
#mine_method_p = "1.957e-01"
mine_method_p = "2.343e-01"
cc=paste("MIC = ",mine_method,"P.value = ",mine_method_p)

xx = normlize(as.numeric(df_method[i,]))
yy = normlize(as.numeric(df_method[j,]))
name = "DeepImpute"

plot(xx, yy, main=name, xlab=i, ylab=j, pch=19) #绘图
legend(x=-0.83,y=1.24, c(cc),bty = "n")

# 绘制第六个图
screen(6)

df_method = df_method7 # 1 2 5 6 7
a = mic_p(df_method[i,],df_method[j,])
mine_method=as.character(a[1])
mine_method_p=as.character(a[2])
mine_method = str_c(substr(mine_method,start=1 ,stop=5),",")
mine_method_p = str_c(substr(mine_method_p,start=1 ,stop=5),"e-",strsplit(mine_method_p,split = "-")[[1]][2])
cc=paste("MIC = ",mine_method,"P.value = ",mine_method_p)

xx = normlize(as.numeric(df_method[i,]))
yy = normlize(as.numeric(df_method[j,]))
name = "scGNGI"

plot(xx, yy, main=name, xlab=i, ylab=j, pch=19) #绘图
legend(x=-0.83,y=1.24, c(cc),bty = "n")

# 结束绘图
dev.off()
# 保存大小 width X height: pdf- 12 8,  eps- 1200 800


#################################################### [3] ####################################################
# 对于不同的细胞基因靶标cell marker，检验其原始值 和 不同方法插补值的相关性(和P值),再画出散点图。
cell_marker <- c("KLF4","NANOG", "SOX2","CD9","CDH11","EFNA2","PRSS50")
########################### plot ########################
# "CD9"   "NANOG" 
# "CD9"  
j = "Raw"  
k = "CD9"

# 归一化
normlize <- function(x){
  a = (x - min(x)) / (max(x) - min(x))
  return(a)
}

cor_p <- function(x_,y_){
  spearman_cor<-c()
  spearman_p<-c()
  x = scale(as.numeric(x_))
  y = scale(as.numeric(y_))
  a=cor(x,y,method = c("spearman"))
  spearman_cor = c(spearman_cor, a)
  
  tmp=cor.test(x,y,method = c("spearman"),alternative = "g")
  #tmp=wilcox.test(x,y)
  spearman_p = c(spearman_p,tmp[["p.value"]]) 
  
  return (list(spearman_cor, spearman_p))
  
}

cor_method<-c()
cor_method_p<-c()

# 分割屏幕
split.screen(c(2, 3))
# 绘制第一个图
screen(1)

library(stringr)
i = "scImpute" #  scImpute SAVER MAGIC ALRA DeepImpute scGNGI 1 2 3 5 6 7
df_method = get(paste0("df_method",as.character(1)))
a = cor_p(df_raw[k,],df_method[k,])
cor_method=as.character(a[1])
cor_method_p=as.character(a[2])

cor_method = str_c(substr(cor_method,start=1 ,stop=5),",")
#cor_method_p = str_c(substr(cor_method_p,start=1 ,stop=5),"e-",strsplit(cor_method_p,split = "-")[[1]][2])
cor_method_p = "0.0e-00"
cc=paste("r = ",cor_method,"P.value = ",cor_method_p)

xx = normlize(as.numeric(df_raw[k,]))
yy = normlize(as.numeric(df_method[k,]))
name = k

plot(xx, yy, main=name, xlab=i, ylab=j, pch=19) #绘图
#legend("topleft", inset=0.001,c(cc),bty = "n",title.adj = 1)
legend(x=-0.6,y=1.24, c(cc),bty = "n")

# 绘制第二个图
screen(2)

i = "SAVER" #  scImpute SAVER MAGIC ALRA DeepImpute scGNGI 1 2 3 5 6 7
df_method = get(paste0("df_method",as.character(2)))
a = cor_p(df_raw[k,],df_method[k,])
cor_method=as.character(a[1])
cor_method_p=as.character(a[2])

cor_method = str_c(substr(cor_method,start=1 ,stop=5),",")
cor_method_p = "0.0e-00"
cc=paste("r = ",cor_method,"P.value = ",cor_method_p)

xx = normlize(as.numeric(df_raw[k,]))
yy = normlize(as.numeric(df_method[k,]))
name = k

plot(xx, yy, main=name, xlab=i, ylab=j, pch=19) #绘图
legend(x=-0.6,y=1.26, c(cc),bty = "n")

# 绘制第三个图
screen(3)

i = "MAGIC" #  scImpute SAVER MAGIC ALRA DeepImpute scGNGI 1 2 3 5 6 7
df_method = get(paste0("df_method",as.character(3)))
a = cor_p(df_raw[k,],df_method[k,])
cor_method=as.character(a[1])
cor_method_p=as.character(a[2])

cor_method = str_c(substr(cor_method,start=1 ,stop=5),",")
cor_method_p = str_c(substr(cor_method_p,start=1 ,stop=5),"e-",strsplit(cor_method_p,split = "-")[[1]][2])
cc=paste("r = ",cor_method,"P.value = ",cor_method_p)

xx = normlize(as.numeric(df_raw[k,]))
yy = normlize(as.numeric(df_method[k,]))
name = k

plot(xx, yy, main=name, xlab=i, ylab=j, pch=19) #绘图
#legend(x=-0.35,y=0.38, c(cc),bty = "n")
cc1=paste("r = ",cor_method)
cc2=paste("P.value = ",cor_method_p)
legend(x=-0.15,y=0.5, c(cc1),bty = "n")
legend(x=-0.15,y=0.38, c(cc2),bty = "n")

# 绘制第四个图
screen(4)

i = "ALRA" #  scImpute SAVER MAGIC ALRA DeepImpute scGNGI 1 2 3 5 6 7
df_method = get(paste0("df_method",as.character(5)))
a = cor_p(df_raw[k,],df_method[k,])
cor_method=as.character(a[1])
cor_method_p=as.character(a[2])

cor_method = str_c(substr(cor_method,start=1 ,stop=5),",")
cor_method_p = str_c(substr(cor_method_p,start=1 ,stop=5),"e-",strsplit(cor_method_p,split = "-")[[1]][2])
cc=paste("r = ",cor_method,"P.value = ",cor_method_p)

xx = normlize(as.numeric(df_raw[k,]))
yy = normlize(as.numeric(df_method[k,]))
name = k

plot(xx, yy, main=name, xlab=i, ylab=j, pch=19) #绘图
legend(x=-0.6,y=1.23, c(cc),bty = "n")


# 绘制第五个图
screen(5)

i = "DeepImpute" #  scImpute SAVER MAGIC ALRA DeepImpute scGNGI 1 2 3 5 6 7
df_method = get(paste0("df_method",as.character(6)))
a = cor_p(df_raw[k,],df_method[k,])
cor_method=as.character(a[1])
cor_method_p=as.character(a[2])

cor_method = str_c(substr(cor_method,start=1 ,stop=5),",")
cor_method_p = str_c(substr(cor_method_p,start=1 ,stop=5),"e-",strsplit(cor_method_p,split = "-")[[1]][2])
cc=paste("r = ",cor_method,"P.value = ",cor_method_p)

xx = normlize(as.numeric(df_raw[k,]))
yy = normlize(as.numeric(df_method[k,]))
name = k

plot(xx, yy, main=name, xlab=i, ylab=j, pch=19) #绘图
legend(x=-0.6,y=1.24, c(cc),bty = "n")


# 绘制第六个图
screen(6)

i = "scGNGI" #  scImpute SAVER MAGIC ALRA DeepImpute scGNGI 1 2 3 5 6 7
df_method = get(paste0("df_method",as.character(7)))
a = cor_p(df_raw[k,],df_method[k,])
cor_method=as.character(a[1])
cor_method_p=as.character(a[2])

cor_method = str_c(substr(cor_method,start=1 ,stop=5),",")
cor_method_p = "0.0e-00"
cc=paste("r = ",cor_method,"P.value = ",cor_method_p)

xx = normlize(as.numeric(df_raw[k,]))
yy = normlize(as.numeric(df_method[k,]))
name = k

plot(xx, yy, main=name, xlab=i, ylab=j, pch=19) #绘图
legend(x=-0.6,y=1.24, c(cc),bty = "n")

# 结束绘图
dev.off()
# 保存大小 width X height: pdf- 12 8,  eps- 1200 800


# "NANOG" 
j = "Raw"  
k = "NANOG"

# 归一化
normlize <- function(x){
  a = (x - min(x)) / (max(x) - min(x))
  return(a)
}

cor_p <- function(x_,y_){
  spearman_cor<-c()
  spearman_p<-c()
  x = scale(as.numeric(x_))
  y = scale(as.numeric(y_))
  a=cor(x,y,method = c("spearman"))
  spearman_cor = c(spearman_cor, a)
  
  tmp=cor.test(x,y,method = c("spearman"),alternative = "g")
  #tmp=wilcox.test(x,y)
  spearman_p = c(spearman_p,tmp[["p.value"]]) 
  
  return (list(spearman_cor, spearman_p))
  
}

cor_method<-c()
cor_method_p<-c()

# 分割屏幕
split.screen(c(2, 3))
# 绘制第一个图
screen(1)

library(stringr)
i = "scImpute" #  scImpute SAVER MAGIC ALRA DeepImpute scGNGI 1 2 3 5 6 7
df_method = get(paste0("df_method",as.character(1)))
a = cor_p(df_raw[k,],df_method[k,])
cor_method=as.character(a[1])
cor_method_p=as.character(a[2])

cor_method = str_c(substr(cor_method,start=1 ,stop=5),",")
#cor_method_p = str_c(substr(cor_method_p,start=1 ,stop=5),"e-",strsplit(cor_method_p,split = "-")[[1]][2])
cor_method_p = "0.0e-00"
cc=paste("r = ",cor_method,"P.value = ",cor_method_p)

xx = normlize(as.numeric(df_raw[k,]))
yy = normlize(as.numeric(df_method[k,]))
name = k

plot(xx, yy, main=name, xlab=i, ylab=j, pch=19) #绘图
#legend("topleft", inset=0.001,c(cc),bty = "n",title.adj = 1)
legend(x=-0.6,y=1.24, c(cc),bty = "n")

# 绘制第二个图
screen(2)

i = "SAVER" #  scImpute SAVER MAGIC ALRA DeepImpute scGNGI 1 2 3 5 6 7
df_method = get(paste0("df_method",as.character(2)))
a = cor_p(df_raw[k,],df_method[k,])
cor_method=as.character(a[1])
cor_method_p=as.character(a[2])

cor_method = str_c(substr(cor_method,start=1 ,stop=5),",")
cor_method_p = "0.0e-00"
cc=paste("r = ",cor_method,"P.value = ",cor_method_p)

xx = normlize(as.numeric(df_raw[k,]))
yy = normlize(as.numeric(df_method[k,]))
name = k

plot(xx, yy, main=name, xlab=i, ylab=j, pch=19) #绘图
legend(x=-0.6,y=1.24, c(cc),bty = "n")

# 绘制第三个图
screen(3)

i = "MAGIC" #  scImpute SAVER MAGIC ALRA DeepImpute scGNGI 1 2 3 5 6 7
df_method = get(paste0("df_method",as.character(3)))
a = cor_p(df_raw[k,],df_method[k,])
cor_method=as.character(a[1])
cor_method_p=as.character(a[2])

cor_method = str_c(substr(cor_method,start=1 ,stop=5),",")
cor_method_p = str_c(substr(cor_method_p,start=1 ,stop=5),"e-",strsplit(cor_method_p,split = "-")[[1]][2])
cc=paste("r = ",cor_method,"P.value = ",cor_method_p)

xx = normlize(as.numeric(df_raw[k,]))
yy = normlize(as.numeric(df_method[k,]))
name = k

plot(xx, yy, main=name, xlab=i, ylab=j, pch=19) #绘图
cc1=paste("r = ",cor_method)
cc2=paste("P.value = ",cor_method_p)
legend(x=-0.16,y=0.5, c(cc1),bty = "n")
legend(x=-0.16,y=0.38, c(cc2),bty = "n")

# 绘制第四个图
screen(4)

i = "ALRA" #  scImpute SAVER MAGIC ALRA DeepImpute scGNGI 1 2 3 5 6 7
df_method = get(paste0("df_method",as.character(5)))
a = cor_p(df_raw[k,],df_method[k,])
cor_method=as.character(a[1])
cor_method_p=as.character(a[2])

cor_method = str_c(substr(cor_method,start=1 ,stop=5),",")
cor_method_p = str_c(substr(cor_method_p,start=1 ,stop=5),"e-",strsplit(cor_method_p,split = "-")[[1]][2])
cc=paste("r = ",cor_method,"P.value = ",cor_method_p)

xx = normlize(as.numeric(df_raw[k,]))
yy = normlize(as.numeric(df_method[k,]))
name = k

plot(xx, yy, main=name, xlab=i, ylab=j, pch=19) #绘图
legend(x=-0.6,y=1.24, c(cc),bty = "n")


# 绘制第五个图
screen(5)

i = "DeepImpute" #  scImpute SAVER MAGIC ALRA DeepImpute scGNGI 1 2 3 5 6 7
df_method = get(paste0("df_method",as.character(6)))
a = cor_p(df_raw[k,],df_method[k,])
cor_method=as.character(a[1])
cor_method_p=as.character(a[2])

cor_method = str_c(substr(cor_method,start=1 ,stop=5),",")
cor_method_p = "0.0e-00"
cc=paste("r = ",cor_method,"P.value = ",cor_method_p)

xx = normlize(as.numeric(df_raw[k,]))
yy = normlize(as.numeric(df_method[k,]))
name = k

plot(xx, yy, main=name, xlab=i, ylab=j, pch=19) #绘图
legend(x=-0.6,y=1.24, c(cc),bty = "n")


# 绘制第六个图
screen(6)

i = "scGNGI" #  scImpute SAVER MAGIC ALRA DeepImpute scGNGI 1 2 3 5 6 7
df_method = get(paste0("df_method",as.character(7)))
a = cor_p(df_raw[k,],df_method[k,])
cor_method=as.character(a[1])
cor_method_p=as.character(a[2])

cor_method = str_c(substr(cor_method,start=1 ,stop=5),",")
cor_method_p = "0.0e-00"
cc=paste("r = ",cor_method,"P.value = ",cor_method_p)

xx = normlize(as.numeric(df_raw[k,]))
yy = normlize(as.numeric(df_method[k,]))
name = k

plot(xx, yy, main=name, xlab=i, ylab=j, pch=19) #绘图
legend(x=-0.6,y=1.24, c(cc),bty = "n")

# 结束绘图
dev.off()
# 保存大小 width X height: pdf- 12 8,  eps- 1200 800
