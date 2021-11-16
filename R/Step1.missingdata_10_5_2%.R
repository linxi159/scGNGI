###########################################
########## creating missing data ##########
###########################################

# 1. GSE75748_sc_cell_type_ec.csv
setwd("~/r_workplace/8.genes_imputation/data")
#gene.expression <- read.csv("./GSE75748/GSE75748_sc_cell_type_ec_proprecessing.csv")
gene.expression <- read.csv("./GSE75748/GSE75748_sc_cell_type_ec_proprecessing_1000_missing82.csv")
rownames(gene.expression) <- gene.expression[,1]
gene.expression <- gene.expression[,-1]

set.seed(1*500+80224523)
gene.expression.missing <- gene.expression
len <- length(gene.expression[gene.expression!=0])
non_zero_rate <- len/ (length(colnames(gene.expression)) * length(rownames(gene.expression)) )
  
missing <- sample(c(0, 1), len, replace=T, prob=c(0.82, 0.18))
gene.expression.missing[gene.expression!=0] <- gene.expression[gene.expression!=0]*missing

write.csv(gene.expression.missing, paste0("./GSE75748/GSE75748_sc_cell_type_ec_proprecessing_1000_missing82", ".csv"))
#write.csv(t(gene.expression.missing), paste0("./GSE75748/GSE75748_sc_cell_type_ec_proprecessing_missing2_t",".csv"))

# 2. GSE90806_RIP-Cre_ARC_GeneCounts.csv
setwd("~/r_workplace/8.genes_imputation/data")
# 去除重复的行
gene.expression <- read.csv("./GSE90806/GSE90806_RIP-Cre_ARC_GeneCounts.csv")
library(dplyr)
gene.expression <- distinct(gene.expression, gene_ID, .keep_all = TRUE)
row.names(gene.expression) <- gene.expression[,1]
gene.expression <- gene.expression[,-1]
write.csv(gene.expression, paste0("./GSE90806/GSE90806_RIP-Cre_ARC_GeneCounts_duplicate_removal", ".csv"))

gene.expression <- read.csv("./GSE90806/GSE90806_RIP-Cre_ARC_GeneCounts_duplicate_removal_proprecessing.csv")
row.names(gene.expression) <- gene.expression[,1]
gene.expression <- gene.expression[,-1]
set.seed(1*500+80224523)
gene.expression.missing <- gene.expression
len <- length(gene.expression[gene.expression!=0])
non_zero_rate <- len / (length(colnames(gene.expression)) * length(rownames(gene.expression)) )

missing <- sample(c(0, 1), len, replace=T, prob=c(0.02, 0.98))
gene.expression.missing[gene.expression!=0] <- gene.expression[gene.expression!=0]*missing

write.csv(gene.expression.missing, paste0("./GSE90806/GSE90806_RIP-Cre_ARC_GeneCounts_duplicate_removal_proprecessing_missing2", ".csv"))
#write.csv(t(gene.expression.missing), paste0("./GSE90806/GSE90806_RIP-Cre_ARC_GeneCounts_duplicate_removal_proprecessing_missing5_t", ".csv"))

# 3. scRNA_data1_group.csv
setwd("~/r_workplace/8.genes_imputation/data")
gene.expression <- read.csv("./scRNA-seq_data_simulation/scRNA_data1.csv")
lbl <- read.csv("./scRNA-seq_data_simulation/scRNA_data1_group.csv")
row.names(gene.expression) <- gene.expression[,1]
gene.expression <- gene.expression[,-1]

# 修改细胞名包含类别并排序
colnames(gene.expression) <- paste(lbl[,2],seq='_',colnames(gene.expression))

gene.expression <- as.data.frame(t(gene.expression))
ls=list(id=row.names(gene.expression))
gene.expression <- cbind(gene.expression,ls) 
gene.expression <- gene.expression[order(gene.expression$id),]
gene.expression <- subset(gene.expression, select=-c(id))
gene.expression <- as.data.frame(t(gene.expression))

write.csv(gene.expression, paste0("./scRNA-seq_data_simulation/scRNA_data1_group1_5", ".csv"))

gene.expression <- read.csv("./scRNA_data1_group1_5/scRNA_data1_group1_5_proprecessing.csv")
row.names(gene.expression) <- gene.expression[,1]
gene.expression <- gene.expression[,-1]
set.seed(1*500+80224523)
gene.expression.missing <- gene.expression
len <- length(gene.expression[gene.expression!=0])
non_zero_rate <- len/ (length(colnames(gene.expression)) * length(rownames(gene.expression)) )

missing <- sample(c(0, 1), len, replace=T, prob=c(0.1, 0.9))
gene.expression.missing[gene.expression!=0] <- gene.expression[gene.expression!=0]*missing

write.csv(gene.expression.missing, paste0("./scRNA_data1_group1_5/scRNA_data1_group1_5_proprecessing_missing10", ".csv"))
#write.csv(t(gene.expression.missing), paste0("./scRNA_data1_group1_5/scRNA_data1_group1_5_proprecessing_missing2_t",".csv"))

# 4. scRNA_data2_group.csv
setwd("~/r_workplace/8.genes_imputation/data")
gene.expression <- read.csv("./scRNA-seq_data_simulation/scRNA_data2.csv")
lbl <- read.csv("./scRNA-seq_data_simulation/scRNA_data2_group.csv")
row.names(gene.expression) <- gene.expression[,1]
gene.expression <- gene.expression[,-1]

# 修改细胞名包含类别并排序
colnames(gene.expression) <- paste(lbl[,2],seq='_',colnames(gene.expression))

gene.expression <- as.data.frame(t(gene.expression))
ls=list(id=row.names(gene.expression))
gene.expression <- cbind(gene.expression,ls) 
gene.expression <- gene.expression[order(gene.expression$id),]
gene.expression <- subset(gene.expression, select=-c(id))
gene.expression <- as.data.frame(t(gene.expression))

write.csv(gene.expression, paste0("./scRNA-seq_data_simulation/scRNA_data2_group1_4", ".csv"))

gene.expression <- read.csv("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing.csv")
row.names(gene.expression) <- gene.expression[,1]
gene.expression <- gene.expression[,-1]
set.seed(1*500+80224523)
gene.expression.missing <- gene.expression
len <- length(gene.expression[gene.expression!=0])
non_zero_rate <- len/ (length(colnames(gene.expression)) * length(rownames(gene.expression)) )

missing <- sample(c(0, 1), len, replace=T, prob=c(0.02, 0.98))
gene.expression.missing[gene.expression!=0] <- gene.expression[gene.expression!=0]*missing

write.csv(gene.expression.missing, paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missing2", ".csv"))
#write.csv(t(gene.expression.missing), paste0("./scRNA_data2_group1_4/scRNA_data2_group1_4_proprecessing_missing10_t",".csv"))

# 5. scRNA_data3_group.csv
setwd("~/r_workplace/8.genes_imputation/data")
gene.expression <- read.csv("./scRNA-seq_data_simulation/scRNA_data3.csv")
lbl <- read.csv("./scRNA-seq_data_simulation/scRNA_data3_group.csv")
row.names(gene.expression) <- gene.expression[,1]
gene.expression <- gene.expression[,-1]

# 修改细胞名包含类别并排序
colnames(gene.expression) <- paste(lbl[,2],seq='_',colnames(gene.expression))

gene.expression <- as.data.frame(t(gene.expression))
ls=list(id=row.names(gene.expression))
gene.expression <- cbind(gene.expression,ls) 
gene.expression <- gene.expression[order(gene.expression$id),]
gene.expression <- subset(gene.expression, select=-c(id))
gene.expression <- as.data.frame(t(gene.expression))

write.csv(gene.expression, paste0("./scRNA-seq_data_simulation/scRNA_data3_group1_3", ".csv"))

gene.expression <- read.csv("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing.csv")
row.names(gene.expression) <- gene.expression[,1]
gene.expression <- gene.expression[,-1]
set.seed(1*500+80224523)
gene.expression.missing <- gene.expression
len <- length(gene.expression[gene.expression!=0])
non_zero_rate <- len/ (length(colnames(gene.expression)) * length(rownames(gene.expression)) )

missing <- sample(c(0, 1), len, replace=T, prob=c(0.35, 0.65))
gene.expression.missing[gene.expression!=0] <- gene.expression[gene.expression!=0]*missing

write.csv(gene.expression.missing, paste0("./scRNA_data3_group1_3/scRNA_data3_group1_3_proprecessing_missing35", ".csv"))





