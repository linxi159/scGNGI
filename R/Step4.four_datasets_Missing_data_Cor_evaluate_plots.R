
####@@@@ missing 10 5 2 @@@@####
rm(list = ls())
#### summary ####
i = 10 # 10 5 2
setwd("~/r_workplace/8.genes_imputation/data")

#### plots ####
library(ggplot2)
library(easyGgplot2)

# correlation_summary
# GSE75748
load(paste0("./GSE75748/results_missing2","/GSE75748_sc_cell_type_ec_proprecessing_missing2","_evaluate_log2.summary"))
df1 <- data.frame(case = rep(1:60, 7), method = c(rep("scImpute", 60), rep("SAVER", 60), rep("MAGIC", 60),rep("DrImpute", 60), rep("ALRA", 60),
                                                  rep("DeepImpute", 60),rep("scGNGI(*)", 60)), value = c(rep(correlation.mean.summary[, 3],60), 
                                                  rep(correlation.mean.summary[, 4],60),rep(correlation.mean.summary[, 5],60),rep(correlation.mean.summary[, 6],60),
                                                  rep(correlation.mean.summary[, 7],60),rep(correlation.mean.summary[, 8],60),rep(correlation.mean.summary[, 9],60)))

load(paste0("./GSE75748/results_missing5","/GSE75748_sc_cell_type_ec_proprecessing_missing5","_evaluate_log2.summary"))
df2 <- data.frame(case = rep(1:60, 7), method = c(rep("scImpute", 60), rep("SAVER", 60), rep("MAGIC", 60),rep("DrImpute", 60), rep("ALRA", 60),
                                                  rep("DeepImpute", 60),rep("scGNGI(*)", 60)), value = c(rep(correlation.mean.summary[, 3],60), 
                                                                                                         rep(correlation.mean.summary[, 4],60),rep(correlation.mean.summary[, 5],60),rep(correlation.mean.summary[, 6],60),
                                                                                                         rep(correlation.mean.summary[, 7],60),rep(correlation.mean.summary[, 8],60),rep(correlation.mean.summary[, 9],60)))

load(paste0("./GSE75748/results_missing10","/GSE75748_sc_cell_type_ec_proprecessing_missing10","_evaluate_log2.summary"))
df3 <- data.frame(case = rep(1:60, 7), method = c(rep("scImpute", 60), rep("SAVER", 60), rep("MAGIC", 60),rep("DrImpute", 60), rep("ALRA", 60),
                                                  rep("DeepImpute", 60),rep("scGNGI(*)", 60)), value = c(rep(correlation.mean.summary[, 3],60), 
                                                                                                         rep(correlation.mean.summary[, 4],60),rep(correlation.mean.summary[, 5],60),rep(correlation.mean.summary[, 6],60),
                                                                                                         rep(correlation.mean.summary[, 7],60),rep(correlation.mean.summary[, 8],60),rep(correlation.mean.summary[, 9],60)))

df <- cbind(missing = c(rep("Missing 2%", 420), rep("Missing 5%", 420), rep("Missing 10%", 420)), rbind(df1, df2, df3))
df$missing <- factor(df$missing, levels = c("Missing 2%", "Missing 5%", "Missing 10%"))
gg1 <- ggplot2.stripchart(data = df, xName = 'method', yName = 'value', groupName = 'method', backgroundColor="white",  
                         groupColors = c('#66B2FF', '#00CC00','#FFAAD4', '#B266FF', '#0000FF','#8000FF', '#FFD4AA'), 
                         xtitle="Imputation Methods", ytitle="Mean Correlation", mainTitle = "Cell Type", removePanelGrid=FALSE, 
                         removePanelBorder=FALSE, boxplotFill="white", showLegend=TRUE, legendTitle = "Method", legendTitleFont = c(10, "bold", "black"), 
                         legendTextFont = c(10, "bold", "black"), setShapeByGroupName=TRUE, addBoxplot=TRUE, faceting=TRUE, facetingVarNames = c("missing"), 
                         facetingDirection="horizontal", mainTitleFont = c(10, "bold", "black"), xtitleFont = c(10, "bold", "black"), size = 4, 
                         ytitleFont = c(10, "bold", "black"), xTickLabelFont = c(10, "bold", "white"), yTickLabelFont = c(10, "bold", "black"), 
                         facetingFont = c(10, "bold", "black"), facetingRect =  list(background = NULL, lineType = NULL, lineColor = NULL, lineSize = NULL)) + scale_shape_manual(values = rep(8:14, len = 7)) 

# GSE90806
load(paste0("./GSE90806/results_missing2","/GSE90806_RIP-Cre_ARC_GeneCounts_duplicate_removal_proprecessing_missing2","_evaluate_log2.summary"))
df1 <- data.frame(case = rep(1:60, 7), method = c(rep("scImpute", 60), rep("SAVER", 60), rep("MAGIC", 60),rep("DrImpute", 60), rep("ALRA", 60),
                                                  rep("DeepImpute", 60),rep("scGNGI(*)", 60)), value = c(rep(correlation.mean.summary[, 3],60), 
                                                                                                         rep(correlation.mean.summary[, 4],60),rep(correlation.mean.summary[, 5],60),rep(correlation.mean.summary[, 6],60),
                                                                                                         rep(correlation.mean.summary[, 7],60),rep(correlation.mean.summary[, 8],60),rep(correlation.mean.summary[, 9],60)))

load(paste0("./GSE90806/results_missing5","/GSE90806_RIP-Cre_ARC_GeneCounts_duplicate_removal_proprecessing_missing5","_evaluate_log2.summary"))
df2 <- data.frame(case = rep(1:60, 7), method = c(rep("scImpute", 60), rep("SAVER", 60), rep("MAGIC", 60),rep("DrImpute", 60), rep("ALRA", 60),
                                                  rep("DeepImpute", 60),rep("scGNGI(*)", 60)), value = c(rep(correlation.mean.summary[, 3],60), 
                                                                                                         rep(correlation.mean.summary[, 4],60),rep(correlation.mean.summary[, 5],60),rep(correlation.mean.summary[, 6],60),
                                                                                                         rep(correlation.mean.summary[, 7],60),rep(correlation.mean.summary[, 8],60),rep(correlation.mean.summary[, 9],60)))

load(paste0("./GSE90806/results_missing10","/GSE90806_RIP-Cre_ARC_GeneCounts_duplicate_removal_proprecessing_missing10","_evaluate_log2.summary"))
df3 <- data.frame(case = rep(1:60, 7), method = c(rep("scImpute", 60), rep("SAVER", 60), rep("MAGIC", 60),rep("DrImpute", 60), rep("ALRA", 60),
                                                  rep("DeepImpute", 60),rep("scGNGI(*)", 60)), value = c(rep(correlation.mean.summary[, 3],60), 
                                                                                                         rep(correlation.mean.summary[, 4],60),rep(correlation.mean.summary[, 5],60),rep(correlation.mean.summary[, 6],60),
                                                                                                         rep(correlation.mean.summary[, 7],60),rep(correlation.mean.summary[, 8],60),rep(correlation.mean.summary[, 9],60)))


df <- cbind(missing = c(rep("Missing 2%", 420), rep("Missing 5%", 420), rep("Missing 10%", 420)), rbind(df1, df2, df3))
df$missing <- factor(df$missing, levels = c("Missing 2%", "Missing 5%", "Missing 10%"))
gg2 <- ggplot2.stripchart(data = df, xName = 'method', yName = 'value', groupName = 'method', backgroundColor="white",  
                          groupColors = c('#66B2FF', '#00CC00','#FFAAD4', '#B266FF', '#0000FF','#8000FF', '#FFD4AA'), 
                          xtitle="Imputation Methods", ytitle="Mean Correlation", mainTitle = "RIP-Cre", removePanelGrid=FALSE, 
                          removePanelBorder=FALSE, boxplotFill="white", showLegend=TRUE, legendTitle = "Method", legendTitleFont = c(10, "bold", "black"), 
                          legendTextFont = c(10, "bold", "black"), setShapeByGroupName=TRUE, addBoxplot=TRUE, faceting=TRUE, facetingVarNames = c("missing"), 
                          facetingDirection="horizontal", mainTitleFont = c(10, "bold", "black"), xtitleFont = c(10, "bold", "black"), size = 4, 
                          ytitleFont = c(10, "bold", "black"), xTickLabelFont = c(10, "bold", "white"), yTickLabelFont = c(10, "bold", "black"), 
                          facetingFont = c(10, "bold", "black"), facetingRect =  list(background = NULL, lineType = NULL, lineColor = NULL, lineSize = NULL)) + scale_shape_manual(values = rep(8:14, len = 7)) 

# scRNA_data1_group1_5
load(paste0("./scRNA_data1_group1_5/results_missing2","/scRNA_data1_group1_5_proprecessing_missing2","_evaluate_log2.summary"))
df1 <- data.frame(case = rep(1:60, 7), method = c(rep("scImpute", 60), rep("SAVER", 60), rep("MAGIC", 60),rep("DrImpute", 60), rep("ALRA", 60),
                                                  rep("DeepImpute", 60),rep("scGNGI(*)", 60)), value = c(rep(correlation.mean.summary[, 3],60), 
                                                                                                         rep(correlation.mean.summary[, 4],60),rep(correlation.mean.summary[, 5],60),rep(correlation.mean.summary[, 6],60),
                                                                                                         rep(correlation.mean.summary[, 7],60),rep(correlation.mean.summary[, 8],60),rep(correlation.mean.summary[, 9],60)))

load(paste0("./scRNA_data1_group1_5/results_missing5","/scRNA_data1_group1_5_proprecessing_missing5","_evaluate_log2.summary"))
df2 <- data.frame(case = rep(1:60, 7), method = c(rep("scImpute", 60), rep("SAVER", 60), rep("MAGIC", 60),rep("DrImpute", 60), rep("ALRA", 60),
                                                  rep("DeepImpute", 60),rep("scGNGI(*)", 60)), value = c(rep(correlation.mean.summary[, 3],60), 
                                                                                                         rep(correlation.mean.summary[, 4],60),rep(correlation.mean.summary[, 5],60),rep(correlation.mean.summary[, 6],60),
                                                                                                         rep(correlation.mean.summary[, 7],60),rep(correlation.mean.summary[, 8],60),rep(correlation.mean.summary[, 9],60)))

load(paste0("./scRNA_data1_group1_5/results_missing10","/scRNA_data1_group1_5_proprecessing_missing10","_evaluate_log2.summary"))
df3 <- data.frame(case = rep(1:60, 7), method = c(rep("scImpute", 60), rep("SAVER", 60), rep("MAGIC", 60),rep("DrImpute", 60), rep("ALRA", 60),
                                                  rep("DeepImpute", 60),rep("scGNGI(*)", 60)), value = c(rep(correlation.mean.summary[, 3],60), 
                                                                                                         rep(correlation.mean.summary[, 4],60),rep(correlation.mean.summary[, 5],60),rep(correlation.mean.summary[, 6],60),
                                                                                                         rep(correlation.mean.summary[, 7],60),rep(correlation.mean.summary[, 8],60),rep(correlation.mean.summary[, 9],60)))


df <- cbind(missing = c(rep("Missing 2%", 420), rep("Missing 5%", 420), rep("Missing 10%", 420)), rbind(df1, df2, df3))
df$missing <- factor(df$missing, levels = c("Missing 2%", "Missing 5%", "Missing 10%"))
gg3 <- ggplot2.stripchart(data = df, xName = 'method', yName = 'value', groupName = 'method', backgroundColor="white",  
                          groupColors = c('#66B2FF', '#00CC00','#FFAAD4', '#B266FF', '#0000FF','#8000FF', '#FFD4AA'), 
                          xtitle="Imputation Methods", ytitle="Mean Correlation", mainTitle = "Simulated Data 1", removePanelGrid=FALSE, 
                          removePanelBorder=FALSE, boxplotFill="white", showLegend=TRUE, legendTitle = "Method", legendTitleFont = c(10, "bold", "black"), 
                          legendTextFont = c(10, "bold", "black"), setShapeByGroupName=TRUE, addBoxplot=TRUE, faceting=TRUE, facetingVarNames = c("missing"), 
                          facetingDirection="horizontal", mainTitleFont = c(10, "bold", "black"), xtitleFont = c(10, "bold", "black"), size = 4, 
                          ytitleFont = c(10, "bold", "black"), xTickLabelFont = c(10, "bold", "white"), yTickLabelFont = c(10, "bold", "black"), 
                          facetingFont = c(10, "bold", "black"), facetingRect =  list(background = NULL, lineType = NULL, lineColor = NULL, lineSize = NULL)) + scale_shape_manual(values = rep(8:14, len = 7)) 


# scRNA_data2_group1_4
load(paste0("./scRNA_data2_group1_4/results_missing2","/scRNA_data2_group1_4_proprecessing_missing2","_evaluate_log2.summary"))
df1 <- data.frame(case = rep(1:60, 7), method = c(rep("scImpute", 60), rep("SAVER", 60), rep("MAGIC", 60),rep("DrImpute", 60), rep("ALRA", 60),
                                                  rep("DeepImpute", 60),rep("scGNGI(*)", 60)), value = c(rep(correlation.mean.summary[, 3],60), 
                                                                                                         rep(correlation.mean.summary[, 4],60),rep(correlation.mean.summary[, 5],60),rep(correlation.mean.summary[, 6],60),
                                                                                                         rep(correlation.mean.summary[, 7],60),rep(correlation.mean.summary[, 8],60),rep(correlation.mean.summary[, 9],60)))

load(paste0("./scRNA_data2_group1_4/results_missing5","/scRNA_data2_group1_4_proprecessing_missing5","_evaluate_log2.summary"))
df2 <- data.frame(case = rep(1:60, 7), method = c(rep("scImpute", 60), rep("SAVER", 60), rep("MAGIC", 60),rep("DrImpute", 60), rep("ALRA", 60),
                                                  rep("DeepImpute", 60),rep("scGNGI(*)", 60)), value = c(rep(correlation.mean.summary[, 3],60), 
                                                                                                         rep(correlation.mean.summary[, 4],60),rep(correlation.mean.summary[, 5],60),rep(correlation.mean.summary[, 6],60),
                                                                                                         rep(correlation.mean.summary[, 7],60),rep(correlation.mean.summary[, 8],60),rep(correlation.mean.summary[, 9],60)))

load(paste0("./scRNA_data2_group1_4/results_missing10","/scRNA_data2_group1_4_proprecessing_missing10","_evaluate_log2.summary"))
df3 <- data.frame(case = rep(1:60, 7), method = c(rep("scImpute", 60), rep("SAVER", 60), rep("MAGIC", 60),rep("DrImpute", 60), rep("ALRA", 60),
                                                  rep("DeepImpute", 60),rep("scGNGI(*)", 60)), value = c(rep(correlation.mean.summary[, 3],60), 
                                                                                                         rep(correlation.mean.summary[, 4],60),rep(correlation.mean.summary[, 5],60),rep(correlation.mean.summary[, 6],60),
                                                                                                         rep(correlation.mean.summary[, 7],60),rep(correlation.mean.summary[, 8],60),rep(correlation.mean.summary[, 9],60)))


df <- cbind(missing = c(rep("Missing 2%", 420), rep("Missing 5%", 420), rep("Missing 10%", 420)), rbind(df1, df2, df3))
df$missing <- factor(df$missing, levels = c("Missing 2%", "Missing 5%", "Missing 10%"))
gg4 <- ggplot2.stripchart(data = df, xName = 'method', yName = 'value', groupName = 'method', backgroundColor="white",  
                          groupColors = c('#66B2FF', '#00CC00','#FFAAD4', '#B266FF', '#0000FF','#8000FF', '#FFD4AA'), 
                          xtitle="Imputation Methods", ytitle="Mean Correlation", mainTitle = "Simulated Data 2", removePanelGrid=FALSE, 
                          removePanelBorder=FALSE, boxplotFill="white", showLegend=TRUE, legendTitle = "Method", legendTitleFont = c(10, "bold", "black"), 
                          legendTextFont = c(10, "bold", "black"), setShapeByGroupName=TRUE, addBoxplot=TRUE, faceting=TRUE, facetingVarNames = c("missing"), 
                          facetingDirection="horizontal", mainTitleFont = c(10, "bold", "black"), xtitleFont = c(10, "bold", "black"), size = 4, 
                          ytitleFont = c(10, "bold", "black"), xTickLabelFont = c(10, "bold", "white"), yTickLabelFont = c(10, "bold", "black"), 
                          facetingFont = c(10, "bold", "black"), facetingRect =  list(background = NULL, lineType = NULL, lineColor = NULL, lineSize = NULL)) + scale_shape_manual(values = rep(8:14, len = 7)) 


ggplot2.multiplot(gg1,gg2,gg3,gg4, cols=1)

#ggsave(ggplot2.multiplot(gg,gg,gg,gg, cols=1), file = paste0("./results/Cor/4datasets_correlation_summary.eps"), width = 12 , height = 16)

# 保存大小 width X height: pdf- 13 11,  eps- 1250 1100(1100 990)

