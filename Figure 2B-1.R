#PART 1:数据清洗################################################################

#载入数据
library(readr)
exp_chemo <- read_csv("exp_chemo.csv")

#制作小鼠的ENSEMBL转SYMBOL表

# 读取GTF文件
library(rtracklayer)
library(dplyr)
gtf <- import('gencode.vM36.annotation.gtf')  # 替换为GTF文件路径

# 转换为数据框并提取需要的信息
gtf_df <- as.data.frame(gtf)

# 提取基因信息
gtf_df <- gtf_df[gtf_df$type == "gene",]

id_conversion <- gtf_df[,c(10,12)]
id_conversion$ENSEMBL <- id_conversion$gene_id
id_conversion$SYMBOL <- id_conversion$gene_name
id_conversion <- id_conversion[,c(3,4)]

# 去除版本号（如果需要）
id_conversion$ENSEMBL <- sub("\\.\\d+$", "", id_conversion$ENSEMBL)

# 更换基因名
exp_chemo <- left_join(exp_chemo, id_conversion, by = "ENSEMBL")

#PART 2: 批量做分析-肺#############################################################
library(dplyr)       
library(tidyr)
library(readr)

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")

# 提取药物-器官对
current_drug_organ <- drug_organ[, c(1, 5)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17), ]
control_organ <- current_drug_organ[14, ]

# 先提取包含 "fei" 的列名并忽略大小写
cols <- grep("fei", names(exp_chemo), value = TRUE, ignore.case = TRUE)

# 对列名进行排序
cols_sorted <- sort(cols)
cols_sorted <- c("ENSEMBL", "SYMBOL", cols_sorted)  # 在前面添加两列

# 根据排序后的列名提取数据
result <- exp_chemo[, cols_sorted]

# 清洗列名中的标识符
colnames(result) <- gsub("\\.fei", "", colnames(result))

result <- result[!is.na(result$SYMBOL), ]
duplicated_rows <- result$SYMBOL[duplicated(result$SYMBOL)]

# 对于 SYMBOL 名相同的情况，保留平均表达量更高的那一行
row_means <- rowMeans(result[, -c(1,2)])
result <- result[order(result$SYMBOL, -row_means), ]
result <- result[!duplicated(result$SYMBOL), ]
rownames(result) <- result$SYMBOL

# 循环批量计算 log2FC，手动计算 log2FC 而非使用 limma 函数
all_results <- list()  # 存储所有比较结果
drug_names <- c()

for (i in 1:length(treat)) {
  # 获取当前处理组和对照组样本名
  treat_list <- treat_organ[i, ]$fei
  treat_list <- unlist(strsplit(treat_list, ", "))
  pbs <- control_organ[control_organ$Drug == control, ]$fei
  pbs <- unlist(strsplit(pbs, ", "))
  
  # 提取表达矩阵（只选择数值列）
  expr_matrix <- result[, c(treat_list, pbs)]
  
  # 计算处理组和对照组的平均表达值
  treat_means <- rowMeans(expr_matrix[, treat_list, drop = FALSE])
  control_means <- rowMeans(expr_matrix[, pbs, drop = FALSE])
  
  # 筛选两组均值大于 0 的基因，避免对数计算中出现无穷或缺失值
  keep_genes <- treat_means > 0 & control_means > 0
  expr_matrix <- expr_matrix[keep_genes, ]
  
  # 确保 expr_matrix 的行名正确
  rownames(expr_matrix) <- rownames(result)[keep_genes]
  
  # 手动计算 log2FC：log2(处理组均值/对照组均值)
  logFC <- log2(treat_means[keep_genes]) - log2(control_means[keep_genes])
  
  # 整理结果数据框，仅保留 logFC 列
  res_df <- data.frame(logFC = logFC)
  rownames(res_df) <- rownames(expr_matrix)
  
  # 保存当前比较的结果
  all_results[[i]] <- res_df
  
  # 获取当前药物名称
  drug_name <- treat_organ[i, ]$Drug
  drug_names <- c(drug_names, drug_name)
}

names(all_results) <- drug_names

# 保存结果，文件命名保持与之前一致
for (i in 1:length(all_results)) {
  drug_name <- names(all_results)[i]
  write.csv(all_results[[i]], 
            file = paste0(drug_name, "_limma_lung_results.csv"),
            quote = FALSE)
}


#PART 3: 批量做分析-血#############################################################
library(dplyr)
library(tidyr)
library(readr)

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")

# 提取药物-器官对 (血)
current_drug_organ <- drug_organ[, c(1,2)]

treat <- current_drug_organ$Drug[c(1:13, 15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13, 15:17), ]
control_organ <- current_drug_organ[14, ]

# 先提取包含 "blood" 的列名并忽略大小写，并增加额外几列
cols <- grep("blood", names(exp_chemo), value = TRUE, ignore.case = TRUE)
cols <- c(cols, "EAP.PBS1", "EAP.PBS2", "EAP.PBS3", "EAP.PBS4")

# 对列名进行排序并在前面添加 "ENSEMBL" 和 "SYMBOL" 列
cols_sorted <- sort(cols)
cols_sorted <- c("ENSEMBL", "SYMBOL", cols_sorted)

# 根据排序后的列名提取数据
result <- exp_chemo[, cols_sorted]

# 清洗列名：去除列名中的 ".blood" 标识符
colnames(result) <- gsub("\\.blood", "", colnames(result))

result <- result[!is.na(result$SYMBOL), ]
# 识别重复的 SYMBOL
duplicated_rows <- result$SYMBOL[duplicated(result$SYMBOL)]

# 对于 SYMBOL 相同，保留平均表达值更高的那一行
row_means <- rowMeans(result[, -c(1,2)])
result <- result[order(result$SYMBOL, -row_means), ]
result <- result[!duplicated(result$SYMBOL), ]
rownames(result) <- result$SYMBOL

# 循环批量手动计算 log2FC
all_results <- list()  # 用于存储所有比较结果
drug_names <- c()

for (i in 1:length(treat)) {
  # 获取当前处理组和对照组样本名（血样数据对应列）
  treat_list <- treat_organ[i, ]$blood
  treat_list <- unlist(strsplit(treat_list, ", "))
  pbs <- control_organ[control_organ$Drug == control, ]$blood
  pbs <- unlist(strsplit(pbs, ", "))
  
  # 提取表达矩阵，选择数值列
  expr_matrix <- result[, c(treat_list, pbs)]
  
  # 计算处理组和对照组的平均表达值
  treat_means <- rowMeans(expr_matrix[, treat_list, drop = FALSE])
  control_means <- rowMeans(expr_matrix[, pbs, drop = FALSE])
  
  # 筛选两组均值均大于 0 的基因，避免对数计算出错
  keep_genes <- treat_means > 0 & control_means > 0
  expr_matrix <- expr_matrix[keep_genes, ]
  
  # 确保 expr_matrix 行名正确
  rownames(expr_matrix) <- rownames(result)[keep_genes]
  
  # 手动计算 log2FC：log2(处理组均值/对照组均值)
  logFC <- log2(treat_means[keep_genes]) - log2(control_means[keep_genes])
  
  # 整理结果，仅保留 logFC 一列
  res_df <- data.frame(logFC = logFC)
  rownames(res_df) <- rownames(expr_matrix)
  
  # 保存当前比较结果
  all_results[[i]] <- res_df
  
  # 保存当前药物名称
  drug_name <- treat_organ[i, ]$Drug
  drug_names <- c(drug_names, drug_name)
}

names(all_results) <- drug_names

# 保存结果文件，命名保持与原来一致
for (i in 1:length(all_results)) {
  drug_name <- names(all_results)[i]
  write.csv(all_results[[i]], 
            file = paste0(drug_name, "_limma_blood_results.csv"),
            quote = FALSE)
}
#PART 4: 批量做分析-心#############################################################
library(dplyr)
library(tidyr)
library(readr)

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")

# 提取药物-器官对 (心)
current_drug_organ <- drug_organ[, c(1,3)]

treat <- current_drug_organ$Drug[c(1:13, 15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13, 15:17), ]
control_organ <- current_drug_organ[14, ]

# 先提取包含 "xin" 的列名并忽略大小写
cols <- grep("xin", names(exp_chemo), value = TRUE, ignore.case = TRUE)

# 对列名进行排序，添加 "ENSEMBL" 和 "SYMBOL" 列
cols_sorted <- sort(cols)
cols_sorted <- c("ENSEMBL", "SYMBOL", cols_sorted)

# 根据排序顺序提取数据
result <- exp_chemo[, cols_sorted]

# 清洗列名：去除列名中的 ".xin" 标识符
colnames(result) <- gsub("\\.xin", "", colnames(result))

result <- result[!is.na(result$SYMBOL), ]
duplicated_rows <- result$SYMBOL[duplicated(result$SYMBOL)]

# 对应 SYMBOL 相同，保留平均表达量更高的那一行
row_means <- rowMeans(result[, -c(1,2)])
result <- result[order(result$SYMBOL, -row_means), ]
result <- result[!duplicated(result$SYMBOL), ]
rownames(result) <- result$SYMBOL

# 循环批量手动计算 log2FC
all_results <- list()  # 用于存储所有比较结果
drug_names <- c()

for (i in 1:length(treat)) {
  # 获取当前处理组和对照组的样本名 (心样数据对应列)
  treat_list <- treat_organ[i, ]$xin
  treat_list <- unlist(strsplit(treat_list, ", "))
  pbs <- control_organ[control_organ$Drug == control, ]$xin
  pbs <- unlist(strsplit(pbs, ", "))
  
  # 提取表达矩阵，只选择数值列
  expr_matrix <- result[, c(treat_list, pbs)]
  
  # 计算处理组和对照组的平均表达值
  treat_means <- rowMeans(expr_matrix[, treat_list, drop = FALSE])
  control_means <- rowMeans(expr_matrix[, pbs, drop = FALSE])
  
  # 筛选两组均值均大于 0 的基因以确保对数计算有效
  keep_genes <- treat_means > 0 & control_means > 0
  expr_matrix <- expr_matrix[keep_genes, ]
  
  # 确保 expr_matrix 行名正确
  rownames(expr_matrix) <- rownames(result)[keep_genes]
  
  # 手动计算 log2FC：log2(处理组均值/对照组均值)
  logFC <- log2(treat_means[keep_genes]) - log2(control_means[keep_genes])
  
  # 整理结果数据框，仅保留 logFC 一列
  res_df <- data.frame(logFC = logFC)
  rownames(res_df) <- rownames(expr_matrix)
  
  # 保存当前比较结果
  all_results[[i]] <- res_df
  
  # 保存当前药物名称
  drug_name <- treat_organ[i, ]$Drug
  drug_names <- c(drug_names, drug_name)
}

names(all_results) <- drug_names

# 保存结果，文件命名保持不变
for (i in 1:length(all_results)) {
  drug_name <- names(all_results)[i]
  write.csv(all_results[[i]], 
            file = paste0(drug_name, "_limma_xin_results.csv"),
            quote = FALSE)
}
#PART 5: 批量做分析-血管#############################################################
library(dplyr)       
library(tidyr)
library(readr)

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")

# 提取药物-器官对
current_drug_organ <- drug_organ[, c(1,4)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17), ]
control_organ <- current_drug_organ[14, ]

# 先提取包含 "guan" 的列名并忽略大小写
cols <- grep("guan", names(exp_chemo), value = TRUE, ignore.case = TRUE)

# 对列名进行排序
cols_sorted <- sort(cols)
cols_sorted <- c("ENSEMBL", "SYMBOL", cols_sorted)  # 在前面添加两列

# 根据排序后的列名提取数据
result <- exp_chemo[, cols_sorted]

# 洗掉列名中的标识符
colnames(result) <- gsub("\\.guan", "", colnames(result))

result <- result[!is.na(result$SYMBOL), ]
duplicated_rows <- result$SYMBOL[duplicated(result$SYMBOL)]

# 对应 SYMBOL 名相同，保留平均表达值更高的那一行
row_means <- rowMeans(result[, -c(1,2)])
result <- result[order(result$SYMBOL, -row_means), ]
result <- result[!duplicated(result$SYMBOL), ]
rownames(result) <- result$SYMBOL

# 循环批量进行差异分析：手动计算 log2FC
all_results <- list()
drug_names <- c()

for (i in 1:length(treat)) {
  # 获取处理组和对照组的样本名
  treat_list <- treat_organ[i, ]$guan
  treat_list <- unlist(strsplit(treat_list, ", "))
  pbs <- control_organ[control_organ$Drug == control, ]$guan
  pbs <- unlist(strsplit(pbs, ", "))
  
  # 提取表达矩阵，只选择数值列
  expr_matrix <- result[, c(treat_list, pbs)]
  
  # 计算处理组和对照组的平均表达值
  treat_means <- rowMeans(expr_matrix[, treat_list, drop = FALSE])
  control_means <- rowMeans(expr_matrix[, pbs, drop = FALSE])
  
  # 筛选两组平均表达值均大于 0 的基因，避免对数计算出错
  keep_genes <- treat_means > 0 & control_means > 0
  expr_matrix <- expr_matrix[keep_genes, ]
  
  # 确保 expr_matrix 继承了 result 的行名
  rownames(expr_matrix) <- rownames(result)[keep_genes]
  
  # 手动计算 log2FC：log2(处理组均值/对照组均值)
  logFC <- log2(treat_means[keep_genes]) - log2(control_means[keep_genes])
  
  # 整理结果数据框，仅保留 logFC 一列
  res_df <- data.frame(logFC = logFC)
  rownames(res_df) <- rownames(expr_matrix)
  
  # 保存当前比较的结果
  all_results[[i]] <- res_df
  
  # 获取当前药物名称
  drug_name <- treat_organ[i, ]$Drug
  drug_names <- c(drug_names, drug_name)
}

names(all_results) <- drug_names

# 保存结果
for (i in 1:length(all_results)) {
  drug_name <- names(all_results)[i]
  write.csv(all_results[[i]], 
            file = paste0(drug_name, "_limma_guan_results.csv"),
            quote = FALSE)
}


#PART 6: 批量做分析-肝#############################################################
library(dplyr)       
library(tidyr)
library(readr)

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")

# 提取药物-器官对
current_drug_organ <- drug_organ[, c(1,6)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17), ]
control_organ <- current_drug_organ[14, ]

# 先提取包含 "gan" 的列名并忽略大小写
cols <- grep("gan", names(exp_chemo), value = TRUE, ignore.case = TRUE)

# 对列名进行排序
cols_sorted <- sort(cols)
cols_sorted <- c("ENSEMBL", "SYMBOL", cols_sorted)  # 在前面添加两列

# 根据排序后的列名提取数据
result <- exp_chemo[, cols_sorted]

# 洗掉列名中的标识符
colnames(result) <- gsub("\\.gan", "", colnames(result))

result <- result[!is.na(result$SYMBOL), ]
duplicated_rows <- result$SYMBOL[duplicated(result$SYMBOL)]

# 对应 SYMBOL 名相同，保留平均表达值更高的那一行
row_means <- rowMeans(result[, -c(1,2)])
result <- result[order(result$SYMBOL, -row_means), ]
result <- result[!duplicated(result$SYMBOL), ]
rownames(result) <- result$SYMBOL

# 循环批量进行差异分析：手动计算 log2FC
all_results <- list()
drug_names <- c()

for (i in 1:length(treat)) {
  # 获取处理组和对照组的样本名
  treat_list <- treat_organ[i, ]$gan
  treat_list <- unlist(strsplit(treat_list, ", "))
  pbs <- control_organ[control_organ$Drug == control, ]$gan
  pbs <- unlist(strsplit(pbs, ", "))
  
  # 提取表达矩阵，只选择数值列
  expr_matrix <- result[, c(treat_list, pbs)]
  
  # 计算处理组和对照组的平均表达值
  treat_means <- rowMeans(expr_matrix[, treat_list, drop = FALSE])
  control_means <- rowMeans(expr_matrix[, pbs, drop = FALSE])
  
  # 筛选两组均值均大于 0 的基因，确保对数计算有效
  keep_genes <- treat_means > 0 & control_means > 0
  expr_matrix <- expr_matrix[keep_genes, ]
  
  # 确保 expr_matrix 继承了 result 的行名
  rownames(expr_matrix) <- rownames(result)[keep_genes]
  
  # 手动计算 log2FC：log2(处理组均值/对照组均值)
  logFC <- log2(treat_means[keep_genes]) - log2(control_means[keep_genes])
  
  # 整理结果，仅保留 logFC 一列
  res_df <- data.frame(logFC = logFC)
  rownames(res_df) <- rownames(expr_matrix)
  
  # 保存当前比较的结果
  all_results[[i]] <- res_df
  
  # 获取当前药物名称
  drug_name <- treat_organ[i, ]$Drug
  drug_names <- c(drug_names, drug_name)
}

names(all_results) <- drug_names

# 保存结果
for (i in 1:length(all_results)) {
  drug_name <- names(all_results)[i]
  write.csv(all_results[[i]], 
            file = paste0(drug_name, "_limma_gan_results.csv"),
            quote = FALSE)
}


#PART 7: 批量做分析-睾丸#############################################################
library(dplyr)       
library(tidyr)
library(readr)

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")

# 提取药物-器官对
current_drug_organ <- drug_organ[, c(1,7)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17), ]
control_organ <- current_drug_organ[14, ]

# 先提取包含 "gao" 的列名并忽略大小写
cols <- grep("gao", names(exp_chemo), value = TRUE, ignore.case = TRUE)

# 对列名进行排序
cols_sorted <- sort(cols)
cols_sorted <- c("ENSEMBL", "SYMBOL", cols_sorted)  # 在前面添加两列

# 根据排序后的列名提取数据
result <- exp_chemo[, cols_sorted]

# 洗掉列名中的标识符
colnames(result) <- gsub("\\.gao", "", colnames(result))

result <- result[!is.na(result$SYMBOL), ]
duplicated_rows <- result$SYMBOL[duplicated(result$SYMBOL)]

# 对应 SYMBOL 名相同，保留平均表达值更高的那一行
row_means <- rowMeans(result[, -c(1,2)])
result <- result[order(result$SYMBOL, -row_means), ]
result <- result[!duplicated(result$SYMBOL), ]
rownames(result) <- result$SYMBOL

# 循环批量进行差异分析：手动计算 log2FC
all_results <- list()
drug_names <- c()

for (i in 1:length(treat)) {
  # 获取处理组和对照组的样本名
  treat_list <- treat_organ[i, ]$gao
  treat_list <- unlist(strsplit(treat_list, ", "))
  pbs <- control_organ[control_organ$Drug == control, ]$gao
  pbs <- unlist(strsplit(pbs, ", "))
  
  # 提取表达矩阵，只选择数值列
  expr_matrix <- result[, c(treat_list, pbs)]
  
  # 计算处理组和对照组的平均表达值
  treat_means <- rowMeans(expr_matrix[, treat_list, drop = FALSE])
  control_means <- rowMeans(expr_matrix[, pbs, drop = FALSE])
  
  # 筛选两组均值均大于 0 的基因，避免对数计算异常
  keep_genes <- treat_means > 0 & control_means > 0
  expr_matrix <- expr_matrix[keep_genes, ]
  
  # 确保 expr_matrix 行名正确
  rownames(expr_matrix) <- rownames(result)[keep_genes]
  
  # 手动计算 log2FC：log2(处理组均值/对照组均值)
  logFC <- log2(treat_means[keep_genes]) - log2(control_means[keep_genes])
  
  # 整理结果数据框，仅保留 logFC 一列
  res_df <- data.frame(logFC = logFC)
  rownames(res_df) <- rownames(expr_matrix)
  
  # 保存当前比较的结果
  all_results[[i]] <- res_df
  
  # 获取当前药物名称
  drug_name <- treat_organ[i, ]$Drug
  drug_names <- c(drug_names, drug_name)
}

names(all_results) <- drug_names

# 保存结果
for (i in 1:length(all_results)) {
  drug_name <- names(all_results)[i]
  write.csv(all_results[[i]], 
            file = paste0(drug_name, "_limma_gao_results.csv"),
            quote = FALSE)
}
#PART 8: 批量做分析-脑#############################################################
library(dplyr)
library(tidyr)
library(readr)

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")

# 提取药物-器官对 (脑)
current_drug_organ <- drug_organ[, c(1,8)]
treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17), ]
control_organ <- current_drug_organ[14, ]

# 先提取包含 "nao" 的列名并忽略大小写
cols <- grep("nao", names(exp_chemo), value = TRUE, ignore.case = TRUE)

# 对列名进行排序
cols_sorted <- sort(cols)
cols_sorted <- c("ENSEMBL", "SYMBOL", cols_sorted)  # 在前面添加两列

# 根据排序后的列名提取数据
result <- exp_chemo[, cols_sorted]

# 洗掉列名中的标识符
colnames(result) <- gsub("\\.nao", "", colnames(result))

result <- result[!is.na(result$SYMBOL), ]
duplicated_rows <- result$SYMBOL[duplicated(result$SYMBOL)]

# 对应 SYMBOL 名相同，保留平均表达量更高的那一行
row_means <- rowMeans(result[, -c(1,2)])
result <- result[order(result$SYMBOL, -row_means), ]
result <- result[!duplicated(result$SYMBOL), ]
rownames(result) <- result$SYMBOL

# 循环批量计算 log2FC
all_results <- list()
drug_names <- c()

for (i in 1:length(treat)) {
  # 获取处理组和对照组的样本名
  treat_list <- treat_organ[i, ]$nao
  treat_list <- unlist(strsplit(treat_list, ", "))
  pbs <- control_organ[control_organ$Drug == control, ]$nao
  pbs <- unlist(strsplit(pbs, ", "))
  
  # 提取表达矩阵（只选择数值列）
  expr_matrix <- result[, c(treat_list, pbs)]
  
  # 计算处理组和对照组的平均表达值
  treat_means <- rowMeans(expr_matrix[, treat_list, drop = FALSE])
  control_means <- rowMeans(expr_matrix[, pbs, drop = FALSE])
  
  # 筛选两组均值大于 0 的基因，确保 log2 计算有效
  keep_genes <- treat_means > 0 & control_means > 0
  expr_matrix <- expr_matrix[keep_genes, ]
  rownames(expr_matrix) <- rownames(result)[keep_genes]
  
  # 手动计算 log2FC：log2(处理组均值) - log2(对照组均值)
  logFC <- log2(treat_means[keep_genes]) - log2(control_means[keep_genes])
  
  # 整理结果数据框，仅保留 logFC 列
  res_df <- data.frame(logFC = logFC)
  rownames(res_df) <- rownames(expr_matrix)
  
  # 保存当前比较的结果
  all_results[[i]] <- res_df
  
  # 获取当前药物名称
  drug_name <- treat_organ[i, ]$Drug
  drug_names <- c(drug_names, drug_name)
}

names(all_results) <- drug_names

# 保存结果文件
for (i in 1:length(all_results)) {
  drug_name <- names(all_results)[i]
  write.csv(all_results[[i]],
            file = paste0(drug_name, "_limma_nao_results.csv"),
            quote = FALSE)
}


#PART 9: 批量做分析-骨#############################################################
library(dplyr)
library(tidyr)
library(readr)

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")

# 提取药物-器官对 (骨)
current_drug_organ <- drug_organ[, c(1,9)]
treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17), ]
control_organ <- current_drug_organ[14, ]

# 先提取包含 "gu" 的列名并忽略大小写
cols <- grep("gu", names(exp_chemo), value = TRUE, ignore.case = TRUE)

# 对列名进行排序
cols_sorted <- sort(cols)
cols_sorted <- c("ENSEMBL", "SYMBOL", cols_sorted)  # 在前面添加两列

# 根据排序后的列名提取数据
result <- exp_chemo[, cols_sorted]

# 洗掉列名中的标识符
colnames(result) <- gsub("\\.gu", "", colnames(result))

result <- result[!is.na(result$SYMBOL), ]
duplicated_rows <- result$SYMBOL[duplicated(result$SYMBOL)]

# 对应 SYMBOL 名相同，保留平均表达量更高的那一行
row_means <- rowMeans(result[, -c(1,2)])
result <- result[order(result$SYMBOL, -row_means), ]
result <- result[!duplicated(result$SYMBOL), ]
rownames(result) <- result$SYMBOL

# 循环批量计算 log2FC
all_results <- list()
drug_names <- c()

for (i in 1:length(treat)) {
  # 获取处理组和对照组的样本名
  treat_list <- treat_organ[i, ]$gu
  treat_list <- unlist(strsplit(treat_list, ", "))
  pbs <- control_organ[control_organ$Drug == control, ]$gu
  pbs <- unlist(strsplit(pbs, ", "))
  
  # 提取表达矩阵，仅选择数值列
  expr_matrix <- result[, c(treat_list, pbs)]
  
  # 计算处理组和对照组的平均表达值
  treat_means <- rowMeans(expr_matrix[, treat_list, drop = FALSE])
  control_means <- rowMeans(expr_matrix[, pbs, drop = FALSE])
  
  # 筛选两组均值大于 0 的基因
  keep_genes <- treat_means > 0 & control_means > 0
  expr_matrix <- expr_matrix[keep_genes, ]
  rownames(expr_matrix) <- rownames(result)[keep_genes]
  
  # 手动计算 log2FC
  logFC <- log2(treat_means[keep_genes]) - log2(control_means[keep_genes])
  
  # 整理为结果数据框（仅 logFC 列）
  res_df <- data.frame(logFC = logFC)
  rownames(res_df) <- rownames(expr_matrix)
  
  all_results[[i]] <- res_df
  
  drug_name <- treat_organ[i, ]$Drug
  drug_names <- c(drug_names, drug_name)
}

names(all_results) <- drug_names

# 保存结果文件
for (i in 1:length(all_results)) {
  drug_name <- names(all_results)[i]
  write.csv(all_results[[i]],
            file = paste0(drug_name, "_limma_gu_results.csv"),
            quote = FALSE)
}


#PART 10: 批量做分析-结肠#############################################################
library(dplyr)
library(tidyr)
library(readr)

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")

# 提取药物-器官对 (结肠)
current_drug_organ <- drug_organ[, c(1,10)]
treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17), ]
control_organ <- current_drug_organ[14, ]

# 先提取包含 "jie" 的列名并忽略大小写
cols <- grep("jie", names(exp_chemo), value = TRUE, ignore.case = TRUE)

# 对列名进行排序
cols_sorted <- sort(cols)
cols_sorted <- c("ENSEMBL", "SYMBOL", cols_sorted)

# 根据排序后的列名提取数据
result <- exp_chemo[, cols_sorted]

# 洗掉列名中的标识符
colnames(result) <- gsub("\\.jie", "", colnames(result))

result <- result[!is.na(result$SYMBOL), ]
duplicated_rows <- result$SYMBOL[duplicated(result$SYMBOL)]

# 保留每个 SYMBOL 平均表达值较高的那一行
row_means <- rowMeans(result[, -c(1,2)])
result <- result[order(result$SYMBOL, -row_means), ]
result <- result[!duplicated(result$SYMBOL), ]
rownames(result) <- result$SYMBOL

# 循环批量计算 log2FC
all_results <- list()
drug_names <- c()

for (i in 1:length(treat)) {
  # 获取处理组和对照组的样本名
  treat_list <- treat_organ[i, ]$jiechang
  treat_list <- unlist(strsplit(treat_list, ", "))
  pbs <- control_organ[control_organ$Drug == control, ]$jiechang
  pbs <- unlist(strsplit(pbs, ", "))
  
  # 提取表达矩阵（只选择数值列）
  expr_matrix <- result[, c(treat_list, pbs)]
  
  # 计算各组平均表达值
  treat_means <- rowMeans(expr_matrix[, treat_list, drop = FALSE])
  control_means <- rowMeans(expr_matrix[, pbs, drop = FALSE])
  
  # 筛选两组均值大于 0 的基因
  keep_genes <- treat_means > 0 & control_means > 0
  expr_matrix <- expr_matrix[keep_genes, ]
  rownames(expr_matrix) <- rownames(result)[keep_genes]
  
  # 计算 log2FC
  logFC <- log2(treat_means[keep_genes]) - log2(control_means[keep_genes])
  
  res_df <- data.frame(logFC = logFC)
  rownames(res_df) <- rownames(expr_matrix)
  
  all_results[[i]] <- res_df
  
  drug_name <- treat_organ[i, ]$Drug
  drug_names <- c(drug_names, drug_name)
}

names(all_results) <- drug_names

# 保存结果文件
for (i in 1:length(all_results)) {
  drug_name <- names(all_results)[i]
  write.csv(all_results[[i]],
            file = paste0(drug_name, "_limma_jiechang_results.csv"),
            quote = FALSE)
}


#PART 11: 批量做分析-肌肉#############################################################
library(dplyr)
library(tidyr)
library(readr)

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")

# 提取药物-器官对 (肌肉)
current_drug_organ <- drug_organ[, c(1,11)]
treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17), ]
control_organ <- current_drug_organ[14, ]

# 先提取包含 "ji" 的列名并忽略大小写
cols <- grep("ji", names(exp_chemo), value = TRUE, ignore.case = TRUE)

# 对列名进行排序
cols_sorted <- sort(cols)
cols_sorted <- c("ENSEMBL", "SYMBOL", cols_sorted)

# 根据排序后的列名提取数据
result <- exp_chemo[, cols_sorted]

# 洗掉列名中的标识符
colnames(result) <- gsub("\\.ji", "", colnames(result))

result <- result[!is.na(result$SYMBOL), ]
duplicated_rows <- result$SYMBOL[duplicated(result$SYMBOL)]

# 保留每个 SYMBOL 中平均表达值较高的行
row_means <- rowMeans(result[, -c(1,2)])
result <- result[order(result$SYMBOL, -row_means), ]
result <- result[!duplicated(result$SYMBOL), ]
rownames(result) <- result$SYMBOL

# 循环批量计算 log2FC
all_results <- list()
drug_names <- c()

for (i in 1:length(treat)) {
  # 获取处理组和对照组的样本名
  treat_list <- treat_organ[i, ]$ji
  treat_list <- unlist(strsplit(treat_list, ", "))
  pbs <- control_organ[control_organ$Drug == control, ]$ji
  pbs <- unlist(strsplit(pbs, ", "))
  
  # 提取表达矩阵（仅选择数值列）
  expr_matrix <- result[, c(treat_list, pbs)]
  
  # 计算两组平均表达值
  treat_means <- rowMeans(expr_matrix[, treat_list, drop = FALSE])
  control_means <- rowMeans(expr_matrix[, pbs, drop = FALSE])
  
  # 筛选两组均值大于 0 的基因
  keep_genes <- treat_means > 0 & control_means > 0
  expr_matrix <- expr_matrix[keep_genes, ]
  rownames(expr_matrix) <- rownames(result)[keep_genes]
  
  # 计算 log2FC
  logFC <- log2(treat_means[keep_genes]) - log2(control_means[keep_genes])
  
  res_df <- data.frame(logFC = logFC)
  rownames(res_df) <- rownames(expr_matrix)
  
  all_results[[i]] <- res_df
  
  drug_name <- treat_organ[i, ]$Drug
  drug_names <- c(drug_names, drug_name)
}

names(all_results) <- drug_names

# 保存结果文件
for (i in 1:length(all_results)) {
  drug_name <- names(all_results)[i]
  write.csv(all_results[[i]],
            file = paste0(drug_name, "_limma_ji_results.csv"),
            quote = FALSE)
}
#PART 12: 批量做分析-胃#############################################################
library(dplyr)
library(tidyr)
library(readr)

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")

# 提取药物-器官对 (胃)
current_drug_organ <- drug_organ[, c(1,12)]
treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17), ]
control_organ <- current_drug_organ[14, ]

# 先提取包含 "wei" 的列名并忽略大小写
cols <- grep("wei", names(exp_chemo), value = TRUE, ignore.case = TRUE)

# 对列名进行排序，并在前面添加"ENSEMBL"和"SYMBOL"
cols_sorted <- sort(cols)
cols_sorted <- c("ENSEMBL", "SYMBOL", cols_sorted)

# 根据排序后的列名提取数据
result <- exp_chemo[, cols_sorted]

# 洗掉列名中的标识符
colnames(result) <- gsub("\\.wei", "", colnames(result))

result <- result[!is.na(result$SYMBOL), ]
# 对于重复的 SYMBOL 保留平均表达值较高的那一行
row_means <- rowMeans(result[, -c(1,2)])
result <- result[order(result$SYMBOL, -row_means), ]
result <- result[!duplicated(result$SYMBOL), ]
rownames(result) <- result$SYMBOL

# 循环批量计算 log2FC（手动计算，不使用 limma 函数）
all_results <- list()
drug_names <- c()

for (i in 1:length(treat)) {
  # 获取处理组和对照组的样本名
  treat_list <- treat_organ[i, ]$wei
  treat_list <- unlist(strsplit(treat_list, ", "))
  pbs <- control_organ[control_organ$Drug == control, ]$wei
  pbs <- unlist(strsplit(pbs, ", "))
  
  # 提取表达矩阵，只选择数值列
  expr_matrix <- result[, c(treat_list, pbs)]
  
  # 计算两组的平均表达值
  treat_means <- rowMeans(expr_matrix[, treat_list, drop = FALSE])
  control_means <- rowMeans(expr_matrix[, pbs, drop = FALSE])
  
  # 筛选两组均值大于0的基因，避免对数计算异常
  keep_genes <- treat_means > 0 & control_means > 0
  expr_matrix <- expr_matrix[keep_genes, ]
  rownames(expr_matrix) <- rownames(result)[keep_genes]
  
  # 手动计算 log2FC：log2(处理组均值) - log2(对照组均值)
  logFC <- log2(treat_means[keep_genes]) - log2(control_means[keep_genes])
  
  res_df <- data.frame(logFC = logFC)
  rownames(res_df) <- rownames(expr_matrix)
  
  all_results[[i]] <- res_df
  
  drug_name <- treat_organ[i, ]$Drug
  drug_names <- c(drug_names, drug_name)
}
names(all_results) <- drug_names

# 保存结果，文件命名保持一致
for (i in 1:length(all_results)) {
  drug_name <- names(all_results)[i]
  write.csv(all_results[[i]],
            file = paste0(drug_name, "_limma_wei_results.csv"),
            quote = FALSE)
}


#PART 13: 批量做分析-皮肤#############################################################
library(dplyr)
library(tidyr)
library(readr)

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")

# 提取药物-器官对 (皮肤)
current_drug_organ <- drug_organ[, c(1,13)]
treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17), ]
control_organ <- current_drug_organ[14, ]

# 先提取包含 "pi" 的列名并忽略大小写
cols <- grep("pi", names(exp_chemo), value = TRUE, ignore.case = TRUE)

# 对列名进行排序，并添加 "ENSEMBL" 与 "SYMBOL"
cols_sorted <- sort(cols)
cols_sorted <- c("ENSEMBL", "SYMBOL", cols_sorted)

# 根据排序后的列名提取数据
result <- exp_chemo[, cols_sorted]

# 洗掉列名中的标识符
colnames(result) <- gsub("\\.pi", "", colnames(result))

result <- result[!is.na(result$SYMBOL), ]
row_means <- rowMeans(result[, -c(1,2)])
result <- result[order(result$SYMBOL, -row_means), ]
result <- result[!duplicated(result$SYMBOL), ]
rownames(result) <- result$SYMBOL

# 循环批量计算 log2FC
all_results <- list()
drug_names <- c()

for (i in 1:length(treat)) {
  # 获取样本名
  treat_list <- treat_organ[i, ]$pi
  treat_list <- unlist(strsplit(treat_list, ", "))
  pbs <- control_organ[control_organ$Drug == control, ]$pi
  pbs <- unlist(strsplit(pbs, ", "))
  
  expr_matrix <- result[, c(treat_list, pbs)]
  treat_means <- rowMeans(expr_matrix[, treat_list, drop = FALSE])
  control_means <- rowMeans(expr_matrix[, pbs, drop = FALSE])
  
  keep_genes <- treat_means > 0 & control_means > 0
  expr_matrix <- expr_matrix[keep_genes, ]
  rownames(expr_matrix) <- rownames(result)[keep_genes]
  
  logFC <- log2(treat_means[keep_genes]) - log2(control_means[keep_genes])
  res_df <- data.frame(logFC = logFC)
  rownames(res_df) <- rownames(expr_matrix)
  
  all_results[[i]] <- res_df
  drug_name <- treat_organ[i, ]$Drug
  drug_names <- c(drug_names, drug_name)
}
names(all_results) <- drug_names

for (i in 1:length(all_results)) {
  drug_name <- names(all_results)[i]
  write.csv(all_results[[i]],
            file = paste0(drug_name, "_limma_pi_results.csv"),
            quote = FALSE)
}


#PART 14: 批量做分析-前列腺#############################################################
library(dplyr)
library(tidyr)
library(readr)

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")

# 提取药物-器官对 (前列腺)
current_drug_organ <- drug_organ[, c(1,14)]
treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17), ]
control_organ <- current_drug_organ[14, ]

# 先提取包含 "qian" 的列名并忽略大小写
cols <- grep("qian", names(exp_chemo), value = TRUE, ignore.case = TRUE)

# 对列名进行排序，并添加 "ENSEMBL" 和 "SYMBOL"
cols_sorted <- sort(cols)
cols_sorted <- c("ENSEMBL", "SYMBOL", cols_sorted)

result <- exp_chemo[, cols_sorted]
colnames(result) <- gsub("\\.qian", "", colnames(result))
result <- result[!is.na(result$SYMBOL), ]
row_means <- rowMeans(result[, -c(1,2)])
result <- result[order(result$SYMBOL, -row_means), ]
result <- result[!duplicated(result$SYMBOL), ]
rownames(result) <- result$SYMBOL

all_results <- list()
drug_names <- c()

for (i in 1:length(treat)) {
  treat_list <- treat_organ[i, ]$qian
  treat_list <- unlist(strsplit(treat_list, ", "))
  pbs <- control_organ[control_organ$Drug == control, ]$qian
  pbs <- unlist(strsplit(pbs, ", "))
  
  expr_matrix <- result[, c(treat_list, pbs)]
  treat_means <- rowMeans(expr_matrix[, treat_list, drop = FALSE])
  control_means <- rowMeans(expr_matrix[, pbs, drop = FALSE])
  
  keep_genes <- treat_means > 0 & control_means > 0
  expr_matrix <- expr_matrix[keep_genes, ]
  rownames(expr_matrix) <- rownames(result)[keep_genes]
  
  logFC <- log2(treat_means[keep_genes]) - log2(control_means[keep_genes])
  res_df <- data.frame(logFC = logFC)
  rownames(res_df) <- rownames(expr_matrix)
  
  all_results[[i]] <- res_df
  drug_name <- treat_organ[i, ]$Drug
  drug_names <- c(drug_names, drug_name)
}
names(all_results) <- drug_names

for (i in 1:length(all_results)) {
  drug_name <- names(all_results)[i]
  write.csv(all_results[[i]],
            file = paste0(drug_name, "_limma_qian_results.csv"),
            quote = FALSE)
}


#PART 15: 批量做分析-肾#############################################################
library(dplyr)
library(tidyr)
library(readr)

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")

# 提取药物-器官对 (肾)
current_drug_organ <- drug_organ[, c(1,15)]
treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17), ]
control_organ <- current_drug_organ[14, ]

# 先提取包含 "shen" 的列名并忽略大小写
cols <- grep("shen", names(exp_chemo), value = TRUE, ignore.case = TRUE)

# 对列名进行排序，并添加 "ENSEMBL" 与 "SYMBOL"
cols_sorted <- sort(cols)
cols_sorted <- c("ENSEMBL", "SYMBOL", cols_sorted)

result <- exp_chemo[, cols_sorted]
colnames(result) <- gsub("\\.shen", "", colnames(result))
result <- result[!is.na(result$SYMBOL), ]
row_means <- rowMeans(result[, -c(1,2)])
result <- result[order(result$SYMBOL, -row_means), ]
result <- result[!duplicated(result$SYMBOL), ]
rownames(result) <- result$SYMBOL

all_results <- list()
drug_names <- c()

for (i in 1:length(treat)) {
  treat_list <- treat_organ[i, ]$shen
  treat_list <- unlist(strsplit(treat_list, ", "))
  pbs <- control_organ[control_organ$Drug == control, ]$shen
  pbs <- unlist(strsplit(pbs, ", "))
  
  expr_matrix <- result[, c(treat_list, pbs)]
  treat_means <- rowMeans(expr_matrix[, treat_list, drop = FALSE])
  control_means <- rowMeans(expr_matrix[, pbs, drop = FALSE])
  
  keep_genes <- treat_means > 0 & control_means > 0
  expr_matrix <- expr_matrix[keep_genes, ]
  rownames(expr_matrix) <- rownames(result)[keep_genes]
  
  logFC <- log2(treat_means[keep_genes]) - log2(control_means[keep_genes])
  res_df <- data.frame(logFC = logFC)
  rownames(res_df) <- rownames(expr_matrix)
  
  all_results[[i]] <- res_df
  drug_name <- treat_organ[i, ]$Drug
  drug_names <- c(drug_names, drug_name)
}
names(all_results) <- drug_names

for (i in 1:length(all_results)) {
  drug_name <- names(all_results)[i]
  write.csv(all_results[[i]],
            file = paste0(drug_name, "_limma_shen_results.csv"),
            quote = FALSE)
}