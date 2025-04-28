#PART 1:数据清洗################################################################

#载入数据
load("E:/9. Chemo_AEs_altas/Analysis/Data/细胞系数据/Chemo_ICI_cell_count_data.Rdata")

shen <- exp_chemo_ici[1]
shen <- as.data.frame(shen)
colnames(shen) <- gsub("X091.1_Count.Rdata.", "", colnames(shen))

nao <- exp_chemo_ici[2]
nao <- as.data.frame(nao)
colnames(nao) <- gsub("X091.10_Count.Rdata.", "", colnames(nao))

ji <- exp_chemo_ici[3]
ji <- as.data.frame(ji)
colnames(ji) <- gsub("X091.11_Count.Rdata.", "", colnames(ji))

gao <- exp_chemo_ici[4]
gao <- as.data.frame(gao)
colnames(gao) <- gsub("X091.12_Count.Rdata.", "", colnames(gao))

fei <- exp_chemo_ici[5]
fei <- as.data.frame(fei)
colnames(fei) <- gsub("X091.13_Count.Rdata.", "", colnames(fei))

wei <- exp_chemo_ici[6]
wei <- as.data.frame(wei)
colnames(wei) <- gsub("X091.2_Count.Rdata.", "", colnames(wei))

gan <- exp_chemo_ici[7]
gan <- as.data.frame(gan)
colnames(gan) <- gsub("X091.3_Count.Rdata.", "", colnames(gan))

gu <- exp_chemo_ici[8]
gu <- as.data.frame(gu)
colnames(gu) <- gsub("X091.4_Count.Rdata.", "", colnames(gu))

xin <- exp_chemo_ici[9]
xin <- as.data.frame(xin)
colnames(xin) <- gsub("X091.5_Count.Rdata.", "", colnames(xin))

guan <- exp_chemo_ici[10]
guan <- as.data.frame(guan)
colnames(guan) <- gsub("X091.6_Count.Rdata.", "", colnames(guan))

jie <- exp_chemo_ici[11]
jie <- as.data.frame(jie)
colnames(jie) <- gsub("X091.7_Count.Rdata.", "", colnames(jie))

qian <- exp_chemo_ici[12]
qian <- as.data.frame(qian)
colnames(qian) <- gsub("X091.8_Count.Rdata.", "", colnames(qian))

pi <- exp_chemo_ici[13]
pi <- as.data.frame(pi)
colnames(pi) <- gsub("X091.9_Count.Rdata.", "", colnames(pi))

#PART 2: 批量做分析-肺#############################################################
library(limma)       # 用于差异分析
library(dplyr)       
library(tidyr)
library(readxl)

# 导入药物-器官分组表
drug_organ <- read_xlsx("drug_organ.xlsx")

# 提取药物-器官对
current_drug_organ <- drug_organ[,c(1,6)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

result <- fei

# 洗掉列名中的标识符
colnames(result) <- gsub("\\_fei", "", colnames(result))

# 循环批量进行差异分析

# 创建一个空的列表来存储所有结果
all_results <- list()
drug_names <- c()

for (i in 1:length(treat)){
  # 获取处理组和对照组的样本名
  treat_list <- treat_organ[i,]$fei
  treat_list <- unlist(strsplit(treat_list, ", "))
  pbs <- control_organ[control_organ$Drug == control,]$fei
  pbs <- unlist(strsplit(pbs, ", "))
  
  # 提取表达矩阵，只选择数值列
  expr_matrix <- result[,c(treat_list,pbs)]
  
  # 计算对照组和处理组的平均值
  treat_means <- rowMeans(expr_matrix[,treat_list,drop=FALSE])
  control_means <- rowMeans(expr_matrix[,pbs,drop=FALSE])
  
  # 筛选两组平均值都不为0的基因
  keep_genes <- treat_means > 0 & control_means > 0
  expr_matrix <- expr_matrix[keep_genes,]
  
  # 确保expr_matrix继承了result的行名
  rownames(expr_matrix) <- rownames(result)[keep_genes]
  
  # 创建设计矩阵
  group <- factor(c(rep("treat", length(treat_list)), 
                    rep("control", length(pbs))))
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  
  # voom转换
  v <- voom(expr_matrix, design)
  
  # 拟合线性模型
  fit <- lmFit(v, design)
  
  # 设置对比
  contrast.matrix <- makeContrasts(treat-control, levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # 提取结果
  res_df <- topTable(fit2, coef=1, number=Inf)
  
  # 添加SYMBOL信息
  res_df$SYMBOL <- rownames(result)[keep_genes][match(rownames(res_df), rownames(expr_matrix))]
  
  # 将SYMBOL设置为行名
  rownames(res_df) <- res_df$SYMBOL
  res_df$SYMBOL <- NULL  # 删除SYMBOL列，因为已经作为行名
  
  # 保存当前比较的结果
  all_results[[i]] <- res_df
  
  # 获取当前药物名称
  drug_name <- treat_organ[i,]$Drug
  drug_names <- c(drug_names, drug_name)
}

names(all_results) <- drug_names

# 保存结果
for (i in 1:length(all_results)) {
  drug_name <- names(all_results)[i]
  write.csv(all_results[[i]], 
            file = paste0(drug_name, "_limma_lung_results.csv"),
            quote = FALSE)
}


#PART 3: 批量做分析-心#############################################################
library(limma)       # 用于差异分析
library(dplyr)       
library(tidyr)
library(readxl)

# 导入药物-器官分组表
drug_organ <- read_xlsx("drug_organ.xlsx")

# 提取药物-器官对
current_drug_organ <- drug_organ[,c(1,10)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

result <- xin

# 洗掉列名中的标识符
colnames(result) <- gsub("\\_xin", "", colnames(result))

# 循环批量进行差异分析

# 创建一个空的列表来存储所有结果
all_results <- list()
drug_names <- c()

for (i in 1:length(treat)){
  # 获取处理组和对照组的样本名
  treat_list <- treat_organ[i,]$xin
  treat_list <- unlist(strsplit(treat_list, ", "))
  pbs <- control_organ[control_organ$Drug == control,]$xin
  pbs <- unlist(strsplit(pbs, ", "))
  
  # 提取表达矩阵，只选择数值列
  expr_matrix <- result[,c(treat_list,pbs)]
  
  # 计算对照组和处理组的平均值
  treat_means <- rowMeans(expr_matrix[,treat_list,drop=FALSE])
  control_means <- rowMeans(expr_matrix[,pbs,drop=FALSE])
  
  # 筛选两组平均值都不为0的基因
  keep_genes <- treat_means > 0 & control_means > 0
  expr_matrix <- expr_matrix[keep_genes,]
  
  # 确保expr_matrix继承了result的行名
  rownames(expr_matrix) <- rownames(result)[keep_genes]
  
  # 创建设计矩阵
  group <- factor(c(rep("treat", length(treat_list)), 
                    rep("control", length(pbs))))
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  
  # voom转换
  v <- voom(expr_matrix, design)
  
  # 拟合线性模型
  fit <- lmFit(v, design)
  
  # 设置对比
  contrast.matrix <- makeContrasts(treat-control, levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # 提取结果
  res_df <- topTable(fit2, coef=1, number=Inf)
  
  # 添加SYMBOL信息
  res_df$SYMBOL <- rownames(result)[keep_genes][match(rownames(res_df), rownames(expr_matrix))]
  
  # 将SYMBOL设置为行名
  rownames(res_df) <- res_df$SYMBOL
  res_df$SYMBOL <- NULL  # 删除SYMBOL列，因为已经作为行名
  
  # 保存当前比较的结果
  all_results[[i]] <- res_df
  
  # 获取当前药物名称
  drug_name <- treat_organ[i,]$Drug
  drug_names <- c(drug_names, drug_name)
}

names(all_results) <- drug_names

# 保存结果
for (i in 1:length(all_results)) {
  drug_name <- names(all_results)[i]
  write.csv(all_results[[i]], 
            file = paste0(drug_name, "_limma_xin_results.csv"),
            quote = FALSE)
}


#PART 4: 批量做分析-肾#############################################################
library(limma)       # 用于差异分析
library(dplyr)       
library(tidyr)
library(readxl)

# 导入药物-器官分组表
drug_organ <- read_xlsx("drug_organ.xlsx")

# 提取药物-器官对
current_drug_organ <- drug_organ[,c(1,2)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

result <- shen

# 洗掉列名中的标识符
colnames(result) <- gsub("\\_shen", "", colnames(result))

# 循环批量进行差异分析

# 创建一个空的列表来存储所有结果
all_results <- list()
drug_names <- c()

for (i in 1:length(treat)){
  # 获取处理组和对照组的样本名
  treat_list <- treat_organ[i,]$shen
  treat_list <- unlist(strsplit(treat_list, ", "))
  pbs <- control_organ[control_organ$Drug == control,]$shen
  pbs <- unlist(strsplit(pbs, ", "))
  
  # 提取表达矩阵，只选择数值列
  expr_matrix <- result[,c(treat_list,pbs)]
  
  # 计算对照组和处理组的平均值
  treat_means <- rowMeans(expr_matrix[,treat_list,drop=FALSE])
  control_means <- rowMeans(expr_matrix[,pbs,drop=FALSE])
  
  # 筛选两组平均值都不为0的基因
  keep_genes <- treat_means > 0 & control_means > 0
  expr_matrix <- expr_matrix[keep_genes,]
  
  # 确保expr_matrix继承了result的行名
  rownames(expr_matrix) <- rownames(result)[keep_genes]
  
  # 创建设计矩阵
  group <- factor(c(rep("treat", length(treat_list)), 
                    rep("control", length(pbs))))
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  
  # voom转换
  v <- voom(expr_matrix, design)
  
  # 拟合线性模型
  fit <- lmFit(v, design)
  
  # 设置对比
  contrast.matrix <- makeContrasts(treat-control, levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # 提取结果
  res_df <- topTable(fit2, coef=1, number=Inf)
  
  # 添加SYMBOL信息
  res_df$SYMBOL <- rownames(result)[keep_genes][match(rownames(res_df), rownames(expr_matrix))]
  
  # 将SYMBOL设置为行名
  rownames(res_df) <- res_df$SYMBOL
  res_df$SYMBOL <- NULL  # 删除SYMBOL列，因为已经作为行名
  
  # 保存当前比较的结果
  all_results[[i]] <- res_df
  
  # 获取当前药物名称
  drug_name <- treat_organ[i,]$Drug
  drug_names <- c(drug_names, drug_name)
}

names(all_results) <- drug_names

# 保存结果
for (i in 1:length(all_results)) {
  drug_name <- names(all_results)[i]
  write.csv(all_results[[i]], 
            file = paste0(drug_name, "_limma_shen_results.csv"),
            quote = FALSE)
}




#PART 5: 批量做分析-脑#############################################################
library(limma)       # 用于差异分析
library(dplyr)       
library(tidyr)
library(readxl)

# 导入药物-器官分组表
drug_organ <- read_xlsx("drug_organ.xlsx")

# 提取药物-器官对
current_drug_organ <- drug_organ[,c(1,3)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

result <- nao

# 洗掉列名中的标识符
colnames(result) <- gsub("\\_jiaozhi", "", colnames(result))

# 循环批量进行差异分析

# 创建一个空的列表来存储所有结果
all_results <- list()
drug_names <- c()

for (i in 1:length(treat)){
  # 获取处理组和对照组的样本名
  treat_list <- treat_organ[i,]$jiaozhi
  treat_list <- unlist(strsplit(treat_list, ", "))
  pbs <- control_organ[control_organ$Drug == control,]$jiaozhi
  pbs <- unlist(strsplit(pbs, ", "))
  
  # 提取表达矩阵，只选择数值列
  expr_matrix <- result[,c(treat_list,pbs)]
  
  # 计算对照组和处理组的平均值
  treat_means <- rowMeans(expr_matrix[,treat_list,drop=FALSE])
  control_means <- rowMeans(expr_matrix[,pbs,drop=FALSE])
  
  # 筛选两组平均值都不为0的基因
  keep_genes <- treat_means > 0 & control_means > 0
  expr_matrix <- expr_matrix[keep_genes,]
  
  # 确保expr_matrix继承了result的行名
  rownames(expr_matrix) <- rownames(result)[keep_genes]
  
  # 创建设计矩阵
  group <- factor(c(rep("treat", length(treat_list)), 
                    rep("control", length(pbs))))
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  
  # voom转换
  v <- voom(expr_matrix, design)
  
  # 拟合线性模型
  fit <- lmFit(v, design)
  
  # 设置对比
  contrast.matrix <- makeContrasts(treat-control, levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # 提取结果
  res_df <- topTable(fit2, coef=1, number=Inf)
  
  # 添加SYMBOL信息
  res_df$SYMBOL <- rownames(result)[keep_genes][match(rownames(res_df), rownames(expr_matrix))]
  
  # 将SYMBOL设置为行名
  rownames(res_df) <- res_df$SYMBOL
  res_df$SYMBOL <- NULL  # 删除SYMBOL列，因为已经作为行名
  
  # 保存当前比较的结果
  all_results[[i]] <- res_df
  
  # 获取当前药物名称
  drug_name <- treat_organ[i,]$Drug
  drug_names <- c(drug_names, drug_name)
}

names(all_results) <- drug_names

# 保存结果
for (i in 1:length(all_results)) {
  drug_name <- names(all_results)[i]
  write.csv(all_results[[i]], 
            file = paste0(drug_name, "_limma_nao_results.csv"),
            quote = FALSE)
}





#PART 6: 批量做分析-肌肉#############################################################
library(limma)       # 用于差异分析
library(dplyr)       
library(tidyr)
library(readxl)

# 导入药物-器官分组表
drug_organ <- read_xlsx("drug_organ.xlsx")

# 提取药物-器官对
current_drug_organ <- drug_organ[,c(1,4)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

result <- ji

# 洗掉列名中的标识符
colnames(result) <- gsub("\\_ji", "", colnames(result))

# 循环批量进行差异分析

# 创建一个空的列表来存储所有结果
all_results <- list()
drug_names <- c()

for (i in 1:length(treat)){
  # 获取处理组和对照组的样本名
  treat_list <- treat_organ[i,]$ji
  treat_list <- unlist(strsplit(treat_list, ", "))
  pbs <- control_organ[control_organ$Drug == control,]$ji
  pbs <- unlist(strsplit(pbs, ", "))
  
  # 提取表达矩阵，只选择数值列
  expr_matrix <- result[,c(treat_list,pbs)]
  
  # 计算对照组和处理组的平均值
  treat_means <- rowMeans(expr_matrix[,treat_list,drop=FALSE])
  control_means <- rowMeans(expr_matrix[,pbs,drop=FALSE])
  
  # 筛选两组平均值都不为0的基因
  keep_genes <- treat_means > 0 & control_means > 0
  expr_matrix <- expr_matrix[keep_genes,]
  
  # 确保expr_matrix继承了result的行名
  rownames(expr_matrix) <- rownames(result)[keep_genes]
  
  # 创建设计矩阵
  group <- factor(c(rep("treat", length(treat_list)), 
                    rep("control", length(pbs))))
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  
  # voom转换
  v <- voom(expr_matrix, design)
  
  # 拟合线性模型
  fit <- lmFit(v, design)
  
  # 设置对比
  contrast.matrix <- makeContrasts(treat-control, levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # 提取结果
  res_df <- topTable(fit2, coef=1, number=Inf)
  
  # 添加SYMBOL信息
  res_df$SYMBOL <- rownames(result)[keep_genes][match(rownames(res_df), rownames(expr_matrix))]
  
  # 将SYMBOL设置为行名
  rownames(res_df) <- res_df$SYMBOL
  res_df$SYMBOL <- NULL  # 删除SYMBOL列，因为已经作为行名
  
  # 保存当前比较的结果
  all_results[[i]] <- res_df
  
  # 获取当前药物名称
  drug_name <- treat_organ[i,]$Drug
  drug_names <- c(drug_names, drug_name)
}

names(all_results) <- drug_names

# 保存结果
for (i in 1:length(all_results)) {
  drug_name <- names(all_results)[i]
  write.csv(all_results[[i]], 
            file = paste0(drug_name, "_limma_ji_results.csv"),
            quote = FALSE)
}





#PART 7: 批量做分析-睾丸#############################################################
library(limma)       # 用于差异分析
library(dplyr)       
library(tidyr)
library(readxl)

# 导入药物-器官分组表
drug_organ <- read_xlsx("drug_organ.xlsx")

# 提取药物-器官对
current_drug_organ <- drug_organ[,c(1,5)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

result <- gao

# 洗掉列名中的标识符
colnames(result) <- gsub("\\_gao", "", colnames(result))

# 循环批量进行差异分析

# 创建一个空的列表来存储所有结果
all_results <- list()
drug_names <- c()

for (i in 1:length(treat)){
  # 获取处理组和对照组的样本名
  treat_list <- treat_organ[i,]$gao
  treat_list <- unlist(strsplit(treat_list, ", "))
  pbs <- control_organ[control_organ$Drug == control,]$gao
  pbs <- unlist(strsplit(pbs, ", "))
  
  # 提取表达矩阵，只选择数值列
  expr_matrix <- result[,c(treat_list,pbs)]
  
  # 计算对照组和处理组的平均值
  treat_means <- rowMeans(expr_matrix[,treat_list,drop=FALSE])
  control_means <- rowMeans(expr_matrix[,pbs,drop=FALSE])
  
  # 筛选两组平均值都不为0的基因
  keep_genes <- treat_means > 0 & control_means > 0
  expr_matrix <- expr_matrix[keep_genes,]
  
  # 确保expr_matrix继承了result的行名
  rownames(expr_matrix) <- rownames(result)[keep_genes]
  
  # 创建设计矩阵
  group <- factor(c(rep("treat", length(treat_list)), 
                    rep("control", length(pbs))))
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  
  # voom转换
  v <- voom(expr_matrix, design)
  
  # 拟合线性模型
  fit <- lmFit(v, design)
  
  # 设置对比
  contrast.matrix <- makeContrasts(treat-control, levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # 提取结果
  res_df <- topTable(fit2, coef=1, number=Inf)
  
  # 添加SYMBOL信息
  res_df$SYMBOL <- rownames(result)[keep_genes][match(rownames(res_df), rownames(expr_matrix))]
  
  # 将SYMBOL设置为行名
  rownames(res_df) <- res_df$SYMBOL
  res_df$SYMBOL <- NULL  # 删除SYMBOL列，因为已经作为行名
  
  # 保存当前比较的结果
  all_results[[i]] <- res_df
  
  # 获取当前药物名称
  drug_name <- treat_organ[i,]$Drug
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





#PART 8: 批量做分析-胃#############################################################
library(limma)       # 用于差异分析
library(dplyr)       
library(tidyr)
library(readxl)

# 导入药物-器官分组表
drug_organ <- read_xlsx("drug_organ.xlsx")

# 提取药物-器官对
current_drug_organ <- drug_organ[,c(1,7)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

result <- wei

# 洗掉列名中的标识符
colnames(result) <- gsub("\\_wei", "", colnames(result))

# 循环批量进行差异分析

# 创建一个空的列表来存储所有结果
all_results <- list()
drug_names <- c()

for (i in 1:length(treat)){
  # 获取处理组和对照组的样本名
  treat_list <- treat_organ[i,]$wei
  treat_list <- unlist(strsplit(treat_list, ", "))
  pbs <- control_organ[control_organ$Drug == control,]$wei
  pbs <- unlist(strsplit(pbs, ", "))
  
  # 提取表达矩阵，只选择数值列
  expr_matrix <- result[,c(treat_list,pbs)]
  
  # 计算对照组和处理组的平均值
  treat_means <- rowMeans(expr_matrix[,treat_list,drop=FALSE])
  control_means <- rowMeans(expr_matrix[,pbs,drop=FALSE])
  
  # 筛选两组平均值都不为0的基因
  keep_genes <- treat_means > 0 & control_means > 0
  expr_matrix <- expr_matrix[keep_genes,]
  
  # 确保expr_matrix继承了result的行名
  rownames(expr_matrix) <- rownames(result)[keep_genes]
  
  # 创建设计矩阵
  group <- factor(c(rep("treat", length(treat_list)), 
                    rep("control", length(pbs))))
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  
  # voom转换
  v <- voom(expr_matrix, design)
  
  # 拟合线性模型
  fit <- lmFit(v, design)
  
  # 设置对比
  contrast.matrix <- makeContrasts(treat-control, levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # 提取结果
  res_df <- topTable(fit2, coef=1, number=Inf)
  
  # 添加SYMBOL信息
  res_df$SYMBOL <- rownames(result)[keep_genes][match(rownames(res_df), rownames(expr_matrix))]
  
  # 将SYMBOL设置为行名
  rownames(res_df) <- res_df$SYMBOL
  res_df$SYMBOL <- NULL  # 删除SYMBOL列，因为已经作为行名
  
  # 保存当前比较的结果
  all_results[[i]] <- res_df
  
  # 获取当前药物名称
  drug_name <- treat_organ[i,]$Drug
  drug_names <- c(drug_names, drug_name)
}

names(all_results) <- drug_names

# 保存结果
for (i in 1:length(all_results)) {
  drug_name <- names(all_results)[i]
  write.csv(all_results[[i]], 
            file = paste0(drug_name, "_limma_wei_results.csv"),
            quote = FALSE)
}





#PART 9: 批量做分析-肝#############################################################
library(limma)       # 用于差异分析
library(dplyr)       
library(tidyr)
library(readxl)

# 导入药物-器官分组表
drug_organ <- read_xlsx("drug_organ.xlsx")

# 提取药物-器官对
current_drug_organ <- drug_organ[,c(1,8)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

result <- gan

# 洗掉列名中的标识符
colnames(result) <- gsub("\\_gan", "", colnames(result))

# 循环批量进行差异分析

# 创建一个空的列表来存储所有结果
all_results <- list()
drug_names <- c()

for (i in 1:length(treat)){
  # 获取处理组和对照组的样本名
  treat_list <- treat_organ[i,]$gan
  treat_list <- unlist(strsplit(treat_list, ", "))
  pbs <- control_organ[control_organ$Drug == control,]$gan
  pbs <- unlist(strsplit(pbs, ", "))
  
  # 提取表达矩阵，只选择数值列
  expr_matrix <- result[,c(treat_list,pbs)]
  
  # 计算对照组和处理组的平均值
  treat_means <- rowMeans(expr_matrix[,treat_list,drop=FALSE])
  control_means <- rowMeans(expr_matrix[,pbs,drop=FALSE])
  
  # 筛选两组平均值都不为0的基因
  keep_genes <- treat_means > 0 & control_means > 0
  expr_matrix <- expr_matrix[keep_genes,]
  
  # 确保expr_matrix继承了result的行名
  rownames(expr_matrix) <- rownames(result)[keep_genes]
  
  # 创建设计矩阵
  group <- factor(c(rep("treat", length(treat_list)), 
                    rep("control", length(pbs))))
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  
  # voom转换
  v <- voom(expr_matrix, design)
  
  # 拟合线性模型
  fit <- lmFit(v, design)
  
  # 设置对比
  contrast.matrix <- makeContrasts(treat-control, levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # 提取结果
  res_df <- topTable(fit2, coef=1, number=Inf)
  
  # 添加SYMBOL信息
  res_df$SYMBOL <- rownames(result)[keep_genes][match(rownames(res_df), rownames(expr_matrix))]
  
  # 将SYMBOL设置为行名
  rownames(res_df) <- res_df$SYMBOL
  res_df$SYMBOL <- NULL  # 删除SYMBOL列，因为已经作为行名
  
  # 保存当前比较的结果
  all_results[[i]] <- res_df
  
  # 获取当前药物名称
  drug_name <- treat_organ[i,]$Drug
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





#PART 10: 批量做分析-骨#############################################################
library(limma)       # 用于差异分析
library(dplyr)       
library(tidyr)
library(readxl)

# 导入药物-器官分组表
drug_organ <- read_xlsx("drug_organ.xlsx")

# 提取药物-器官对
current_drug_organ <- drug_organ[,c(1,9)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

result <- gu

# 洗掉列名中的标识符
colnames(result) <- gsub("\\_gu", "", colnames(result))

# 循环批量进行差异分析

# 创建一个空的列表来存储所有结果
all_results <- list()
drug_names <- c()

for (i in 1:length(treat)){
  # 获取处理组和对照组的样本名
  treat_list <- treat_organ[i,]$gu
  treat_list <- unlist(strsplit(treat_list, ", "))
  pbs <- control_organ[control_organ$Drug == control,]$gu
  pbs <- unlist(strsplit(pbs, ", "))
  
  # 提取表达矩阵，只选择数值列
  expr_matrix <- result[,c(treat_list,pbs)]
  
  # 计算对照组和处理组的平均值
  treat_means <- rowMeans(expr_matrix[,treat_list,drop=FALSE])
  control_means <- rowMeans(expr_matrix[,pbs,drop=FALSE])
  
  # 筛选两组平均值都不为0的基因
  keep_genes <- treat_means > 0 & control_means > 0
  expr_matrix <- expr_matrix[keep_genes,]
  
  # 确保expr_matrix继承了result的行名
  rownames(expr_matrix) <- rownames(result)[keep_genes]
  
  # 创建设计矩阵
  group <- factor(c(rep("treat", length(treat_list)), 
                    rep("control", length(pbs))))
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  
  # voom转换
  v <- voom(expr_matrix, design)
  
  # 拟合线性模型
  fit <- lmFit(v, design)
  
  # 设置对比
  contrast.matrix <- makeContrasts(treat-control, levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # 提取结果
  res_df <- topTable(fit2, coef=1, number=Inf)
  
  # 添加SYMBOL信息
  res_df$SYMBOL <- rownames(result)[keep_genes][match(rownames(res_df), rownames(expr_matrix))]
  
  # 将SYMBOL设置为行名
  rownames(res_df) <- res_df$SYMBOL
  res_df$SYMBOL <- NULL  # 删除SYMBOL列，因为已经作为行名
  
  # 保存当前比较的结果
  all_results[[i]] <- res_df
  
  # 获取当前药物名称
  drug_name <- treat_organ[i,]$Drug
  drug_names <- c(drug_names, drug_name)
}

names(all_results) <- drug_names

# 保存结果
for (i in 1:length(all_results)) {
  drug_name <- names(all_results)[i]
  write.csv(all_results[[i]], 
            file = paste0(drug_name, "_limma_gu_results.csv"),
            quote = FALSE)
}





#PART 11: 批量做分析-血管#############################################################
library(limma)       # 用于差异分析
library(dplyr)       
library(tidyr)
library(readxl)

# 导入药物-器官分组表
drug_organ <- read_xlsx("drug_organ.xlsx")

# 提取药物-器官对
current_drug_organ <- drug_organ[,c(1,11)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

result <- guan

# 洗掉列名中的标识符
colnames(result) <- gsub("\\_neipi", "", colnames(result))

# 循环批量进行差异分析

# 创建一个空的列表来存储所有结果
all_results <- list()
drug_names <- c()

for (i in 1:length(treat)){
  # 获取处理组和对照组的样本名
  treat_list <- treat_organ[i,]$neipi
  treat_list <- unlist(strsplit(treat_list, ", "))
  pbs <- control_organ[control_organ$Drug == control,]$neipi
  pbs <- unlist(strsplit(pbs, ", "))
  
  # 提取表达矩阵，只选择数值列
  expr_matrix <- result[,c(treat_list,pbs)]
  
  # 计算对照组和处理组的平均值
  treat_means <- rowMeans(expr_matrix[,treat_list,drop=FALSE])
  control_means <- rowMeans(expr_matrix[,pbs,drop=FALSE])
  
  # 筛选两组平均值都不为0的基因
  keep_genes <- treat_means > 0 & control_means > 0
  expr_matrix <- expr_matrix[keep_genes,]
  
  # 确保expr_matrix继承了result的行名
  rownames(expr_matrix) <- rownames(result)[keep_genes]
  
  # 创建设计矩阵
  group <- factor(c(rep("treat", length(treat_list)), 
                    rep("control", length(pbs))))
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  
  # voom转换
  v <- voom(expr_matrix, design)
  
  # 拟合线性模型
  fit <- lmFit(v, design)
  
  # 设置对比
  contrast.matrix <- makeContrasts(treat-control, levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # 提取结果
  res_df <- topTable(fit2, coef=1, number=Inf)
  
  # 添加SYMBOL信息
  res_df$SYMBOL <- rownames(result)[keep_genes][match(rownames(res_df), rownames(expr_matrix))]
  
  # 将SYMBOL设置为行名
  rownames(res_df) <- res_df$SYMBOL
  res_df$SYMBOL <- NULL  # 删除SYMBOL列，因为已经作为行名
  
  # 保存当前比较的结果
  all_results[[i]] <- res_df
  
  # 获取当前药物名称
  drug_name <- treat_organ[i,]$Drug
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





#PART 12: 批量做分析-结肠#############################################################
library(limma)       # 用于差异分析
library(dplyr)       
library(tidyr)
library(readxl)

# 导入药物-器官分组表
drug_organ <- read_xlsx("drug_organ.xlsx")

# 提取药物-器官对
current_drug_organ <- drug_organ[,c(1,12)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

result <- jie

# 洗掉列名中的标识符
colnames(result) <- gsub("\\_jie", "", colnames(result))

# 循环批量进行差异分析

# 创建一个空的列表来存储所有结果
all_results <- list()
drug_names <- c()

for (i in 1:length(treat)){
  # 获取处理组和对照组的样本名
  treat_list <- treat_organ[i,]$jie
  treat_list <- unlist(strsplit(treat_list, ", "))
  pbs <- control_organ[control_organ$Drug == control,]$jie
  pbs <- unlist(strsplit(pbs, ", "))
  
  # 提取表达矩阵，只选择数值列
  expr_matrix <- result[,c(treat_list,pbs)]
  
  # 计算对照组和处理组的平均值
  treat_means <- rowMeans(expr_matrix[,treat_list,drop=FALSE])
  control_means <- rowMeans(expr_matrix[,pbs,drop=FALSE])
  
  # 筛选两组平均值都不为0的基因
  keep_genes <- treat_means > 0 & control_means > 0
  expr_matrix <- expr_matrix[keep_genes,]
  
  # 确保expr_matrix继承了result的行名
  rownames(expr_matrix) <- rownames(result)[keep_genes]
  
  # 创建设计矩阵
  group <- factor(c(rep("treat", length(treat_list)), 
                    rep("control", length(pbs))))
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  
  # voom转换
  v <- voom(expr_matrix, design)
  
  # 拟合线性模型
  fit <- lmFit(v, design)
  
  # 设置对比
  contrast.matrix <- makeContrasts(treat-control, levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # 提取结果
  res_df <- topTable(fit2, coef=1, number=Inf)
  
  # 添加SYMBOL信息
  res_df$SYMBOL <- rownames(result)[keep_genes][match(rownames(res_df), rownames(expr_matrix))]
  
  # 将SYMBOL设置为行名
  rownames(res_df) <- res_df$SYMBOL
  res_df$SYMBOL <- NULL  # 删除SYMBOL列，因为已经作为行名
  
  # 保存当前比较的结果
  all_results[[i]] <- res_df
  
  # 获取当前药物名称
  drug_name <- treat_organ[i,]$Drug
  drug_names <- c(drug_names, drug_name)
}

names(all_results) <- drug_names

# 保存结果
for (i in 1:length(all_results)) {
  drug_name <- names(all_results)[i]
  write.csv(all_results[[i]], 
            file = paste0(drug_name, "_limma_jie_results.csv"),
            quote = FALSE)
}





#PART 13: 批量做分析-前列腺#############################################################
library(limma)       # 用于差异分析
library(dplyr)       
library(tidyr)
library(readxl)

# 导入药物-器官分组表
drug_organ <- read_xlsx("drug_organ.xlsx")

# 提取药物-器官对
current_drug_organ <- drug_organ[,c(1,13)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

result <- qian

# 洗掉列名中的标识符
colnames(result) <- gsub("\\_qian", "", colnames(result))

# 循环批量进行差异分析

# 创建一个空的列表来存储所有结果
all_results <- list()
drug_names <- c()

for (i in 1:length(treat)){
  # 获取处理组和对照组的样本名
  treat_list <- treat_organ[i,]$qian
  treat_list <- unlist(strsplit(treat_list, ", "))
  pbs <- control_organ[control_organ$Drug == control,]$qian
  pbs <- unlist(strsplit(pbs, ", "))
  
  # 提取表达矩阵，只选择数值列
  expr_matrix <- result[,c(treat_list,pbs)]
  
  # 计算对照组和处理组的平均值
  treat_means <- rowMeans(expr_matrix[,treat_list,drop=FALSE])
  control_means <- rowMeans(expr_matrix[,pbs,drop=FALSE])
  
  # 筛选两组平均值都不为0的基因
  keep_genes <- treat_means > 0 & control_means > 0
  expr_matrix <- expr_matrix[keep_genes,]
  
  # 确保expr_matrix继承了result的行名
  rownames(expr_matrix) <- rownames(result)[keep_genes]
  
  # 创建设计矩阵
  group <- factor(c(rep("treat", length(treat_list)), 
                    rep("control", length(pbs))))
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  
  # voom转换
  v <- voom(expr_matrix, design)
  
  # 拟合线性模型
  fit <- lmFit(v, design)
  
  # 设置对比
  contrast.matrix <- makeContrasts(treat-control, levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # 提取结果
  res_df <- topTable(fit2, coef=1, number=Inf)
  
  # 添加SYMBOL信息
  res_df$SYMBOL <- rownames(result)[keep_genes][match(rownames(res_df), rownames(expr_matrix))]
  
  # 将SYMBOL设置为行名
  rownames(res_df) <- res_df$SYMBOL
  res_df$SYMBOL <- NULL  # 删除SYMBOL列，因为已经作为行名
  
  # 保存当前比较的结果
  all_results[[i]] <- res_df
  
  # 获取当前药物名称
  drug_name <- treat_organ[i,]$Drug
  drug_names <- c(drug_names, drug_name)
}

names(all_results) <- drug_names

# 保存结果
for (i in 1:length(all_results)) {
  drug_name <- names(all_results)[i]
  write.csv(all_results[[i]], 
            file = paste0(drug_name, "_limma_qian_results.csv"),
            quote = FALSE)
}





#PART 14: 批量做分析-皮肤#############################################################
library(limma)       # 用于差异分析
library(dplyr)       
library(tidyr)
library(readxl)

# 导入药物-器官分组表
drug_organ <- read_xlsx("drug_organ.xlsx")

# 提取药物-器官对
current_drug_organ <- drug_organ[,c(1,14)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

result <- pi

# 洗掉列名中的标识符
colnames(result) <- gsub("\\_biao", "", colnames(result))

# 循环批量进行差异分析

# 创建一个空的列表来存储所有结果
all_results <- list()
drug_names <- c()

for (i in 1:length(treat)){
  # 获取处理组和对照组的样本名
  treat_list <- treat_organ[i,]$biao
  treat_list <- unlist(strsplit(treat_list, ", "))
  pbs <- control_organ[control_organ$Drug == control,]$biao
  pbs <- unlist(strsplit(pbs, ", "))
  
  # 提取表达矩阵，只选择数值列
  expr_matrix <- result[,c(treat_list,pbs)]
  
  # 计算对照组和处理组的平均值
  treat_means <- rowMeans(expr_matrix[,treat_list,drop=FALSE])
  control_means <- rowMeans(expr_matrix[,pbs,drop=FALSE])
  
  # 筛选两组平均值都不为0的基因
  keep_genes <- treat_means > 0 & control_means > 0
  expr_matrix <- expr_matrix[keep_genes,]
  
  # 确保expr_matrix继承了result的行名
  rownames(expr_matrix) <- rownames(result)[keep_genes]
  
  # 创建设计矩阵
  group <- factor(c(rep("treat", length(treat_list)), 
                    rep("control", length(pbs))))
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  
  # voom转换
  v <- voom(expr_matrix, design)
  
  # 拟合线性模型
  fit <- lmFit(v, design)
  
  # 设置对比
  contrast.matrix <- makeContrasts(treat-control, levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # 提取结果
  res_df <- topTable(fit2, coef=1, number=Inf)
  
  # 添加SYMBOL信息
  res_df$SYMBOL <- rownames(result)[keep_genes][match(rownames(res_df), rownames(expr_matrix))]
  
  # 将SYMBOL设置为行名
  rownames(res_df) <- res_df$SYMBOL
  res_df$SYMBOL <- NULL  # 删除SYMBOL列，因为已经作为行名
  
  # 保存当前比较的结果
  all_results[[i]] <- res_df
  
  # 获取当前药物名称
  drug_name <- treat_organ[i,]$Drug
  drug_names <- c(drug_names, drug_name)
}

names(all_results) <- drug_names

# 保存结果
for (i in 1:length(all_results)) {
  drug_name <- names(all_results)[i]
  write.csv(all_results[[i]], 
            file = paste0(drug_name, "_limma_pi_results.csv"),
            quote = FALSE)
}




