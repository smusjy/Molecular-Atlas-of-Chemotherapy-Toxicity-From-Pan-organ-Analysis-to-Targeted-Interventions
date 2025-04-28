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

# 首先记录每个数据框的原始行数
original_rows <- c(
  shen = nrow(shen),
  nao = nrow(nao),
  ji = nrow(ji),
  gao = nrow(gao),
  fei = nrow(fei),
  wei = nrow(wei),
  gan = nrow(gan),
  gu = nrow(gu),
  xin = nrow(xin),
  guan = nrow(guan),
  jie = nrow(jie),
  qian = nrow(qian),
  pi = nrow(pi)
)

# 添加基因名称列
shen$gene <- rownames(shen)
nao$gene <- rownames(nao)
ji$gene <- rownames(ji)
gao$gene <- rownames(gao)
fei$gene <- rownames(fei)
wei$gene <- rownames(wei)
gan$gene <- rownames(gan)
gu$gene <- rownames(gu)
xin$gene <- rownames(xin)
guan$gene <- rownames(guan)
jie$gene <- rownames(jie)
qian$gene <- rownames(qian)
pi$gene <- rownames(pi)

# 合并数据
merged_data <- Reduce(function(x, y) merge(x, y, by="gene", all=TRUE), 
                      list(shen, nao, ji, gao, fei, wei, gan, 
                           gu, xin, guan, jie, qian, pi))

# 设置行名并删除gene列
rownames(merged_data) <- merged_data$gene
merged_data$gene <- NULL

# 验证步骤
print("原始数据框行数：")
print(original_rows)
print("合并后数据框行数：")
print(nrow(merged_data))

#PART 2: 批量做分析-肺#############################################################
library(limma)       # 用于差异分析
library(dplyr)       
library(tidyr)
library(readxl)

# 导入药物-器官分组表
drug_organ <- read_xlsx("drug_organ.xlsx")

# 提取药物-器官对
current_drug_organ <- drug_organ[,c(1,6)]

# 获取所有处理组和对照组
treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

result <- fei

# 洗掉列名中的标识符
colnames(result) <- gsub("\\_fei", "", colnames(result))

# 获取所有处理组样本和对照组样本
all_treat_samples <- unlist(strsplit(as.character(treat_organ$fei), ", "))
control_samples <- unlist(strsplit(as.character(control_organ$fei), ", "))

# 提取表达矩阵
expr_matrix <- result[,c(all_treat_samples, control_samples)]

# 计算处理组和对照组的平均值
treat_means <- rowMeans(expr_matrix[,all_treat_samples,drop=FALSE])
control_means <- rowMeans(expr_matrix[,control_samples,drop=FALSE])

# 筛选两组平均值都不为0的基因
keep_genes <- treat_means > 0 & control_means > 0
expr_matrix <- expr_matrix[keep_genes,]

# 确保expr_matrix继承了result的行名
rownames(expr_matrix) <- rownames(result)[keep_genes]

# 创建设计矩阵
group <- factor(c(rep("treat", length(all_treat_samples)), 
                  rep("control", length(control_samples))))
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

# 手动计算log2FC
treat_means_filtered <- treat_means[keep_genes]
control_means_filtered <- control_means[keep_genes]
manual_log2FC <- log2(treat_means_filtered) - log2(control_means_filtered)

# 将手动计算的log2FC添加到结果中
res_df$manual_log2FC <- manual_log2FC[match(rownames(res_df), names(manual_log2FC))]

# 添加SYMBOL信息
res_df$SYMBOL <- rownames(result)[keep_genes][match(rownames(res_df), rownames(expr_matrix))]

# 将SYMBOL设置为行名
rownames(res_df) <- res_df$SYMBOL
res_df$SYMBOL <- NULL  # 删除SYMBOL列，因为已经作为行名

# 保存结果
write.csv(res_df, file = "all_drugs_vs_PBS_limma_lung_results.csv", quote = FALSE)
#PART 3: 批量做分析-肾#############################################################
library(limma)       # 用于差异分析
library(dplyr)       
library(tidyr)
library(readxl)

# 导入药物-器官分组表
drug_organ <- read_xlsx("drug_organ.xlsx")

# 提取药物-器官对
current_drug_organ <- drug_organ[,c(1,2)]

# 获取所有处理组和对照组
treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

result <- shen

# 洗掉列名中的标识符
colnames(result) <- gsub("\\_shen", "", colnames(result))

# 获取所有处理组样本和对照组样本
all_treat_samples <- unlist(strsplit(as.character(treat_organ$shen), ", "))
control_samples <- unlist(strsplit(as.character(control_organ$shen), ", "))

# 提取表达矩阵
expr_matrix <- result[,c(all_treat_samples, control_samples)]

# 计算处理组和对照组的平均值
treat_means <- rowMeans(expr_matrix[,all_treat_samples,drop=FALSE])
control_means <- rowMeans(expr_matrix[,control_samples,drop=FALSE])

# 筛选两组平均值都不为0的基因
keep_genes <- treat_means > 0 & control_means > 0
expr_matrix <- expr_matrix[keep_genes,]

# 确保expr_matrix继承了result的行名
rownames(expr_matrix) <- rownames(result)[keep_genes]

# 手动计算log2FC
manual_log2FC <- log2(treat_means[keep_genes]) - log2(control_means[keep_genes])

# 创建设计矩阵
group <- factor(c(rep("treat", length(all_treat_samples)), 
                  rep("control", length(control_samples))))
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
res_df$manual_log2FC <- manual_log2FC[match(rownames(res_df), names(manual_log2FC))]

# 添加SYMBOL信息
res_df$SYMBOL <- rownames(result)[keep_genes][match(rownames(res_df), rownames(expr_matrix))]

# 将SYMBOL设置为行名
rownames(res_df) <- res_df$SYMBOL
res_df$SYMBOL <- NULL  # 删除SYMBOL列，因为已经作为行名

# 保存结果
write.csv(res_df, file = "all_drugs_vs_PBS_limma_shen_results.csv", quote = FALSE)

#PART 4: 批量做分析-脑#############################################################
library(limma)       # 用于差异分析
library(dplyr)       
library(tidyr)
library(readxl)

# 导入药物-器官分组表
drug_organ <- read_xlsx("drug_organ.xlsx")

# 提取药物-器官对
current_drug_organ <- drug_organ[,c(1,3)]

# 获取所有处理组和对照组
treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

result <- nao

# 洗掉列名中的标识符
colnames(result) <- gsub("\\_jiaozhi", "", colnames(result))

# 获取所有处理组样本和对照组样本
all_treat_samples <- unlist(strsplit(as.character(treat_organ$jiaozhi), ", "))
control_samples <- unlist(strsplit(as.character(control_organ$jiaozhi), ", "))

# 提取表达矩阵
expr_matrix <- result[,c(all_treat_samples, control_samples)]

# 计算处理组和对照组的平均值
treat_means <- rowMeans(expr_matrix[,all_treat_samples,drop=FALSE])
control_means <- rowMeans(expr_matrix[,control_samples,drop=FALSE])

# 筛选两组平均值都不为0的基因
keep_genes <- treat_means > 0 & control_means > 0
expr_matrix <- expr_matrix[keep_genes,]

# 确保expr_matrix继承了result的行名
rownames(expr_matrix) <- rownames(result)[keep_genes]

# 手动计算log2FC
manual_log2FC <- log2(treat_means[keep_genes]) - log2(control_means[keep_genes])

# 创建设计矩阵
group <- factor(c(rep("treat", length(all_treat_samples)), 
                  rep("control", length(control_samples))))
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
res_df$manual_log2FC <- manual_log2FC[match(rownames(res_df), names(manual_log2FC))]

# 添加SYMBOL信息
res_df$SYMBOL <- rownames(result)[keep_genes][match(rownames(res_df), rownames(expr_matrix))]

# 将SYMBOL设置为行名
rownames(res_df) <- res_df$SYMBOL
res_df$SYMBOL <- NULL  # 删除SYMBOL列，因为已经作为行名

# 保存结果
write.csv(res_df, file = "all_drugs_vs_PBS_limma_nao_results.csv", quote = FALSE)

#PART 5: 批量做分析-肌肉#############################################################
library(limma)       # 用于差异分析
library(dplyr)       
library(tidyr)
library(readxl)

# 导入药物-器官分组表
drug_organ <- read_xlsx("drug_organ.xlsx")

# 提取药物-器官对
current_drug_organ <- drug_organ[,c(1,4)]

# 获取所有处理组和对照组
treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

result <- ji

# 洗掉列名中的标识符
colnames(result) <- gsub("\\_ji", "", colnames(result))

# 获取所有处理组样本和对照组样本
all_treat_samples <- unlist(strsplit(as.character(treat_organ$ji), ", "))
control_samples <- unlist(strsplit(as.character(control_organ$ji), ", "))

# 提取表达矩阵
expr_matrix <- result[,c(all_treat_samples, control_samples)]

# 计算处理组和对照组的平均值
treat_means <- rowMeans(expr_matrix[,all_treat_samples,drop=FALSE])
control_means <- rowMeans(expr_matrix[,control_samples,drop=FALSE])

# 筛选两组平均值都不为0的基因
keep_genes <- treat_means > 0 & control_means > 0
expr_matrix <- expr_matrix[keep_genes,]

# 确保expr_matrix继承了result的行名
rownames(expr_matrix) <- rownames(result)[keep_genes]

# 手动计算log2FC
manual_log2FC <- log2(treat_means[keep_genes]) - log2(control_means[keep_genes])

# 创建设计矩阵
group <- factor(c(rep("treat", length(all_treat_samples)), 
                  rep("control", length(control_samples))))
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
res_df$manual_log2FC <- manual_log2FC[match(rownames(res_df), names(manual_log2FC))]

# 添加SYMBOL信息
res_df$SYMBOL <- rownames(result)[keep_genes][match(rownames(res_df), rownames(expr_matrix))]

# 将SYMBOL设置为行名
rownames(res_df) <- res_df$SYMBOL
res_df$SYMBOL <- NULL  # 删除SYMBOL列，因为已经作为行名

# 保存结果
write.csv(res_df, file = "all_drugs_vs_PBS_limma_ji_results.csv", quote = FALSE)
#PART 6: 批量做分析-睾丸#############################################################
library(limma)       # 用于差异分析
library(dplyr)       
library(tidyr)
library(readxl)

# 导入药物-器官分组表
drug_organ <- read_xlsx("drug_organ.xlsx")

# 提取药物-器官对
current_drug_organ <- drug_organ[,c(1,5)]

# 获取所有处理组和对照组
treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

result <- gao

# 洗掉列名中的标识符
colnames(result) <- gsub("\\_gao", "", colnames(result))

# 获取所有处理组样本和对照组样本
all_treat_samples <- unlist(strsplit(as.character(treat_organ$gao), ", "))
control_samples <- unlist(strsplit(as.character(control_organ$gao), ", "))

# 提取表达矩阵
expr_matrix <- result[,c(all_treat_samples, control_samples)]

# 计算处理组和对照组的平均值
treat_means <- rowMeans(expr_matrix[,all_treat_samples,drop=FALSE])
control_means <- rowMeans(expr_matrix[,control_samples,drop=FALSE])

# 筛选两组平均值都不为0的基因
keep_genes <- treat_means > 0 & control_means > 0
expr_matrix <- expr_matrix[keep_genes,]

# 确保expr_matrix继承了result的行名
rownames(expr_matrix) <- rownames(result)[keep_genes]

# 手动计算log2FC
manual_log2FC <- log2(treat_means[keep_genes]) - log2(control_means[keep_genes])

# 创建设计矩阵
group <- factor(c(rep("treat", length(all_treat_samples)), 
                  rep("control", length(control_samples))))
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
res_df$manual_log2FC <- manual_log2FC[match(rownames(res_df), names(manual_log2FC))]

# 添加SYMBOL信息
res_df$SYMBOL <- rownames(result)[keep_genes][match(rownames(res_df), rownames(expr_matrix))]

# 将SYMBOL设置为行名
rownames(res_df) <- res_df$SYMBOL
res_df$SYMBOL <- NULL  # 删除SYMBOL列，因为已经作为行名

# 保存结果
write.csv(res_df, file = "all_drugs_vs_PBS_limma_gao_results.csv", quote = FALSE)

#PART 7: 批量做分析-胃#############################################################
library(limma)       # 用于差异分析
library(dplyr)       
library(tidyr)
library(readxl)

# 导入药物-器官分组表
drug_organ <- read_xlsx("drug_organ.xlsx")

# 提取药物-器官对
current_drug_organ <- drug_organ[,c(1,7)]

# 获取所有处理组和对照组
treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

result <- wei

# 洗掉列名中的标识符
colnames(result) <- gsub("\\_wei", "", colnames(result))

# 获取所有处理组样本和对照组样本
all_treat_samples <- unlist(strsplit(as.character(treat_organ$wei), ", "))
control_samples <- unlist(strsplit(as.character(control_organ$wei), ", "))

# 提取表达矩阵
expr_matrix <- result[,c(all_treat_samples, control_samples)]

# 计算处理组和对照组的平均值
treat_means <- rowMeans(expr_matrix[,all_treat_samples,drop=FALSE])
control_means <- rowMeans(expr_matrix[,control_samples,drop=FALSE])

# 筛选两组平均值都不为0的基因
keep_genes <- treat_means > 0 & control_means > 0
expr_matrix <- expr_matrix[keep_genes,]

# 确保expr_matrix继承了result的行名
rownames(expr_matrix) <- rownames(result)[keep_genes]

# 手动计算log2FC
manual_log2FC <- log2(treat_means[keep_genes]) - log2(control_means[keep_genes])

# 创建设计矩阵
group <- factor(c(rep("treat", length(all_treat_samples)), 
                  rep("control", length(control_samples))))
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
res_df$manual_log2FC <- manual_log2FC[match(rownames(res_df), names(manual_log2FC))]

# 添加SYMBOL信息
res_df$SYMBOL <- rownames(result)[keep_genes][match(rownames(res_df), rownames(expr_matrix))]

# 将SYMBOL设置为行名
rownames(res_df) <- res_df$SYMBOL
res_df$SYMBOL <- NULL  # 删除SYMBOL列，因为已经作为行名

# 保存结果
write.csv(res_df, file = "all_drugs_vs_PBS_limma_wei_results.csv", quote = FALSE)

#PART 8: 批量做分析-肝#############################################################
library(limma)       # 用于差异分析
library(dplyr)       
library(tidyr)
library(readxl)

# 导入药物-器官分组表
drug_organ <- read_xlsx("drug_organ.xlsx")

# 提取药物-器官对
current_drug_organ <- drug_organ[,c(1,8)]

# 获取所有处理组和对照组
treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

result <- gan

# 洗掉列名中的标识符
colnames(result) <- gsub("\\_gan", "", colnames(result))

# 获取所有处理组样本和对照组样本
all_treat_samples <- unlist(strsplit(as.character(treat_organ$gan), ", "))
control_samples <- unlist(strsplit(as.character(control_organ$gan), ", "))

# 提取表达矩阵
expr_matrix <- result[,c(all_treat_samples, control_samples)]

# 计算处理组和对照组的平均值
treat_means <- rowMeans(expr_matrix[,all_treat_samples,drop=FALSE])
control_means <- rowMeans(expr_matrix[,control_samples,drop=FALSE])

# 筛选两组平均值都不为0的基因
keep_genes <- treat_means > 0 & control_means > 0
expr_matrix <- expr_matrix[keep_genes,]

# 确保expr_matrix继承了result的行名
rownames(expr_matrix) <- rownames(result)[keep_genes]

# 手动计算log2FC
manual_log2FC <- log2(treat_means[keep_genes]) - log2(control_means[keep_genes])

# 创建设计矩阵
group <- factor(c(rep("treat", length(all_treat_samples)), 
                  rep("control", length(control_samples))))
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
res_df$manual_log2FC <- manual_log2FC[match(rownames(res_df), names(manual_log2FC))]

# 添加SYMBOL信息
res_df$SYMBOL <- rownames(result)[keep_genes][match(rownames(res_df), rownames(expr_matrix))]

# 将SYMBOL设置为行名
rownames(res_df) <- res_df$SYMBOL
res_df$SYMBOL <- NULL  # 删除SYMBOL列，因为已经作为行名

# 保存结果
write.csv(res_df, file = "all_drugs_vs_PBS_limma_gan_results.csv", quote = FALSE)
#PART 9: 批量做分析-骨#############################################################
library(limma)       # 用于差异分析
library(dplyr)       
library(tidyr)
library(readxl)

# 导入药物-器官分组表
drug_organ <- read_xlsx("drug_organ.xlsx")

# 提取药物-器官对
current_drug_organ <- drug_organ[,c(1,9)]

# 获取所有处理组和对照组
treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

result <- gu

# 洗掉列名中的标识符
colnames(result) <- gsub("\\_gu", "", colnames(result))

# 获取所有处理组样本和对照组样本
all_treat_samples <- unlist(strsplit(as.character(treat_organ$gu), ", "))
control_samples <- unlist(strsplit(as.character(control_organ$gu), ", "))

# 提取表达矩阵
expr_matrix <- result[,c(all_treat_samples, control_samples)]

# 计算处理组和对照组的平均值
treat_means <- rowMeans(expr_matrix[,all_treat_samples,drop=FALSE])
control_means <- rowMeans(expr_matrix[,control_samples,drop=FALSE])

# 筛选两组平均值都不为0的基因
keep_genes <- treat_means > 0 & control_means > 0
expr_matrix <- expr_matrix[keep_genes,]

# 确保expr_matrix继承了result的行名
rownames(expr_matrix) <- rownames(result)[keep_genes]

# 手动计算log2FC
manual_log2FC <- log2(treat_means[keep_genes]) - log2(control_means[keep_genes])

# 创建设计矩阵
group <- factor(c(rep("treat", length(all_treat_samples)), 
                  rep("control", length(control_samples))))
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
res_df$manual_log2FC <- manual_log2FC[match(rownames(res_df), names(manual_log2FC))]

# 添加SYMBOL信息
res_df$SYMBOL <- rownames(result)[keep_genes][match(rownames(res_df), rownames(expr_matrix))]

# 将SYMBOL设置为行名
rownames(res_df) <- res_df$SYMBOL
res_df$SYMBOL <- NULL  # 删除SYMBOL列，因为已经作为行名

# 保存结果
write.csv(res_df, file = "all_drugs_vs_PBS_limma_gu_results.csv", quote = FALSE)

#PART 10: 批量做分析-心#############################################################
library(limma)       # 用于差异分析
library(dplyr)       
library(tidyr)
library(readxl)

# 导入药物-器官分组表
drug_organ <- read_xlsx("drug_organ.xlsx")

# 提取药物-器官对
current_drug_organ <- drug_organ[,c(1,10)]

# 获取所有处理组和对照组
treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

result <- xin

# 洗掉列名中的标识符
colnames(result) <- gsub("\\_xin", "", colnames(result))

# 获取所有处理组样本和对照组样本
all_treat_samples <- unlist(strsplit(as.character(treat_organ$xin), ", "))
control_samples <- unlist(strsplit(as.character(control_organ$xin), ", "))

# 提取表达矩阵
expr_matrix <- result[,c(all_treat_samples, control_samples)]

# 计算处理组和对照组的平均值
treat_means <- rowMeans(expr_matrix[,all_treat_samples,drop=FALSE])
control_means <- rowMeans(expr_matrix[,control_samples,drop=FALSE])

# 筛选两组平均值都不为0的基因
keep_genes <- treat_means > 0 & control_means > 0
expr_matrix <- expr_matrix[keep_genes,]

# 确保expr_matrix继承了result的行名
rownames(expr_matrix) <- rownames(result)[keep_genes]

# 手动计算log2FC
manual_log2FC <- log2(treat_means[keep_genes]) - log2(control_means[keep_genes])

# 创建设计矩阵
group <- factor(c(rep("treat", length(all_treat_samples)), 
                  rep("control", length(control_samples))))
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
res_df$manual_log2FC <- manual_log2FC[match(rownames(res_df), names(manual_log2FC))]

# 添加SYMBOL信息
res_df$SYMBOL <- rownames(result)[keep_genes][match(rownames(res_df), rownames(expr_matrix))]

# 将SYMBOL设置为行名
rownames(res_df) <- res_df$SYMBOL
res_df$SYMBOL <- NULL  # 删除SYMBOL列，因为已经作为行名

# 保存结果
write.csv(res_df, file = "all_drugs_vs_PBS_limma_xin_results.csv", quote = FALSE)

#PART 11: 批量做分析-血管#############################################################
library(limma)       # 用于差异分析
library(dplyr)       
library(tidyr)
library(readxl)

# 导入药物-器官分组表
drug_organ <- read_xlsx("drug_organ.xlsx")

# 提取药物-器官对
current_drug_organ <- drug_organ[,c(1,11)]

# 获取所有处理组和对照组
treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

result <- guan

# 洗掉列名中的标识符
colnames(result) <- gsub("\\_neipi", "", colnames(result))

# 获取所有处理组样本和对照组样本
all_treat_samples <- unlist(strsplit(as.character(treat_organ$neipi), ", "))
control_samples <- unlist(strsplit(as.character(control_organ$neipi), ", "))

# 提取表达矩阵
expr_matrix <- result[,c(all_treat_samples, control_samples)]

# 计算处理组和对照组的平均值
treat_means <- rowMeans(expr_matrix[,all_treat_samples,drop=FALSE])
control_means <- rowMeans(expr_matrix[,control_samples,drop=FALSE])

# 筛选两组平均值都不为0的基因
keep_genes <- treat_means > 0 & control_means > 0
expr_matrix <- expr_matrix[keep_genes,]

# 确保expr_matrix继承了result的行名
rownames(expr_matrix) <- rownames(result)[keep_genes]

# 手动计算log2FC
manual_log2FC <- log2(treat_means[keep_genes]) - log2(control_means[keep_genes])

# 创建设计矩阵
group <- factor(c(rep("treat", length(all_treat_samples)), 
                  rep("control", length(control_samples))))
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
res_df$manual_log2FC <- manual_log2FC[match(rownames(res_df), names(manual_log2FC))]

# 添加SYMBOL信息
res_df$SYMBOL <- rownames(result)[keep_genes][match(rownames(res_df), rownames(expr_matrix))]

# 将SYMBOL设置为行名
rownames(res_df) <- res_df$SYMBOL
res_df$SYMBOL <- NULL  # 删除SYMBOL列，因为已经作为行名

# 保存结果
write.csv(res_df, file = "all_drugs_vs_PBS_limma_guan_results.csv", quote = FALSE)
#PART 12: 批量做分析-结肠#############################################################
library(limma)       # 用于差异分析
library(dplyr)       
library(tidyr)
library(readxl)

# 导入药物-器官分组表
drug_organ <- read_xlsx("drug_organ.xlsx")

# 提取药物-器官对
current_drug_organ <- drug_organ[,c(1,12)]

# 获取所有处理组和对照组
treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

result <- jie

# 洗掉列名中的标识符
colnames(result) <- gsub("\\_jie", "", colnames(result))

# 获取所有处理组样本和对照组样本
all_treat_samples <- unlist(strsplit(as.character(treat_organ$jie), ", "))
control_samples <- unlist(strsplit(as.character(control_organ$jie), ", "))

# 提取表达矩阵
expr_matrix <- result[,c(all_treat_samples, control_samples)]

# 计算处理组和对照组的平均值
treat_means <- rowMeans(expr_matrix[,all_treat_samples,drop=FALSE])
control_means <- rowMeans(expr_matrix[,control_samples,drop=FALSE])

# 筛选两组平均值都不为0的基因
keep_genes <- treat_means > 0 & control_means > 0
expr_matrix <- expr_matrix[keep_genes,]

# 确保expr_matrix继承了result的行名
rownames(expr_matrix) <- rownames(result)[keep_genes]

# 手动计算log2FC
manual_log2FC <- log2(treat_means[keep_genes]) - log2(control_means[keep_genes])

# 创建设计矩阵
group <- factor(c(rep("treat", length(all_treat_samples)), 
                  rep("control", length(control_samples))))
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
res_df$manual_log2FC <- manual_log2FC[match(rownames(res_df), names(manual_log2FC))]

# 添加SYMBOL信息
res_df$SYMBOL <- rownames(result)[keep_genes][match(rownames(res_df), rownames(expr_matrix))]

# 将SYMBOL设置为行名
rownames(res_df) <- res_df$SYMBOL
res_df$SYMBOL <- NULL  # 删除SYMBOL列，因为已经作为行名

# 保存结果
write.csv(res_df, file = "all_drugs_vs_PBS_limma_jie_results.csv", quote = FALSE)

#PART 13: 批量做分析-前列腺#############################################################
library(limma)       # 用于差异分析
library(dplyr)       
library(tidyr)
library(readxl)

# 导入药物-器官分组表
drug_organ <- read_xlsx("drug_organ.xlsx")

# 提取药物-器官对
current_drug_organ <- drug_organ[,c(1,13)]

# 获取所有处理组和对照组
treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

result <- qian

# 洗掉列名中的标识符
colnames(result) <- gsub("\\_qian", "", colnames(result))

# 获取所有处理组样本和对照组样本
all_treat_samples <- unlist(strsplit(as.character(treat_organ$qian), ", "))
control_samples <- unlist(strsplit(as.character(control_organ$qian), ", "))

# 提取表达矩阵
expr_matrix <- result[,c(all_treat_samples, control_samples)]

# 计算处理组和对照组的平均值
treat_means <- rowMeans(expr_matrix[,all_treat_samples,drop=FALSE])
control_means <- rowMeans(expr_matrix[,control_samples,drop=FALSE])

# 筛选两组平均值都不为0的基因
keep_genes <- treat_means > 0 & control_means > 0
expr_matrix <- expr_matrix[keep_genes,]

# 确保expr_matrix继承了result的行名
rownames(expr_matrix) <- rownames(result)[keep_genes]

# 手动计算log2FC
manual_log2FC <- log2(treat_means[keep_genes]) - log2(control_means[keep_genes])

# 创建设计矩阵
group <- factor(c(rep("treat", length(all_treat_samples)), 
                  rep("control", length(control_samples))))
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
res_df$manual_log2FC <- manual_log2FC[match(rownames(res_df), names(manual_log2FC))]

# 添加SYMBOL信息
res_df$SYMBOL <- rownames(result)[keep_genes][match(rownames(res_df), rownames(expr_matrix))]

# 将SYMBOL设置为行名
rownames(res_df) <- res_df$SYMBOL
res_df$SYMBOL <- NULL  # 删除SYMBOL列，因为已经作为行名

# 保存结果
write.csv(res_df, file = "all_drugs_vs_PBS_limma_qian_results.csv", quote = FALSE)

#PART 14: 批量做分析-皮肤#############################################################
library(limma)       # 用于差异分析
library(dplyr)       
library(tidyr)
library(readxl)

# 导入药物-器官分组表
drug_organ <- read_xlsx("drug_organ.xlsx")

# 提取药物-器官对
current_drug_organ <- drug_organ[,c(1,14)]

# 获取所有处理组和对照组
treat <- current_drug_organ$Drug[c(1:13,15:17)]
control <- current_drug_organ$Drug[14]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

result <- pi

# 洗掉列名中的标识符
colnames(result) <- gsub("\\_biao", "", colnames(result))

# 获取所有处理组样本和对照组样本
all_treat_samples <- unlist(strsplit(as.character(treat_organ$biao), ", "))
control_samples <- unlist(strsplit(as.character(control_organ$biao), ", "))

# 提取表达矩阵
expr_matrix <- result[,c(all_treat_samples, control_samples)]

# 计算处理组和对照组的平均值
treat_means <- rowMeans(expr_matrix[,all_treat_samples,drop=FALSE])
control_means <- rowMeans(expr_matrix[,control_samples,drop=FALSE])

# 筛选两组平均值都不为0的基因
keep_genes <- treat_means > 0 & control_means > 0
expr_matrix <- expr_matrix[keep_genes,]

# 确保expr_matrix继承了result的行名
rownames(expr_matrix) <- rownames(result)[keep_genes]

# 手动计算log2FC
manual_log2FC <- log2(treat_means[keep_genes]) - log2(control_means[keep_genes])

# 创建设计矩阵
group <- factor(c(rep("treat", length(all_treat_samples)), 
                  rep("control", length(control_samples))))
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
res_df$manual_log2FC <- manual_log2FC[match(rownames(res_df), names(manual_log2FC))]

# 添加SYMBOL信息
res_df$SYMBOL <- rownames(result)[keep_genes][match(rownames(res_df), rownames(expr_matrix))]

# 将SYMBOL设置为行名
rownames(res_df) <- res_df$SYMBOL
res_df$SYMBOL <- NULL  # 删除SYMBOL列，因为已经作为行名

# 保存结果
write.csv(res_df, file = "all_drugs_vs_PBS_limma_pi_results.csv", quote = FALSE)