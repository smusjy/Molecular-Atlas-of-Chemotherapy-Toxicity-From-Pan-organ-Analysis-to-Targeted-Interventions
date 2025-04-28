# PART 1: 数据清洗 ###############################################################
library(readr)
library(dplyr)

# 载入表达数据
exp_chemo <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/Protein_coding/exp_chemo.csv")

# 制作小鼠的 ENSEMBL 转 SYMBOL 表
library(rtracklayer)
gtf <- import('E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/Protein_coding/gencode.vM36.annotation.gtf')  # 替换为GTF文件路径

# 转换为数据框并提取需要的信息，只保留基因层级的记录
gtf_df <- as.data.frame(gtf)
gtf_df <- gtf_df[gtf_df$type == "gene",]

# 假设第10、11、12列为 gene_id, gene_name, gene_type
id_conversion <- gtf_df[, c(10, 11, 12)]
id_conversion$ENSEMBL <- id_conversion$gene_id
id_conversion$SYMBOL  <- id_conversion$gene_name
id_conversion <- id_conversion[, c(2, 4, 5)]
# 仅保留 protein_coding 基因
id_conversion <- id_conversion[id_conversion$gene_type == "protein_coding", ]
id_conversion <- id_conversion[, c(2, 3)]
# 去除版本号（如果需要）
id_conversion$ENSEMBL <- sub("\\.\\d+$", "", id_conversion$ENSEMBL)

# 用 SYMBOL 更换 exp_chemo 中的 ENSEMBL id，并去除未匹配上 SYMBOL 的行
exp_chemo <- left_join(exp_chemo, id_conversion, by = "ENSEMBL")
exp_chemo <- exp_chemo[!is.na(exp_chemo$SYMBOL),]

# 合并重复的基因（假设表达数值在列2到列1359）
merged_data <- exp_chemo
numeric_cols <- merged_data[, c(2:1359)]
gene_names  <- merged_data$SYMBOL
sums <- rowsum(numeric_cols, gene_names)
counts <- as.vector(table(gene_names)[rownames(sums)])
merged_vals <- sums / counts
merged_data <- data.frame(SYMBOL = rownames(merged_vals), merged_vals)


# PART 2: 按器官计算泛化疗处理组与PBS对照的 log2FC ###############################
library(tidyr)

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")

# 现在把治疗组固定为所有非PBS行，PBS行作为对照
treat_data <- drug_organ %>% filter(Drug != "PBS")
control_data <- drug_organ %>% filter(Drug == "PBS")

# 假设药物-器官表中第一列是 Drug，后续各列代表不同器官
organ_cols <- names(drug_organ)[-1]

# 对每个器官进行循环
for(organ in organ_cols) {
  cat("开始分析器官:", organ, "\n")
  
  # 根据器官名称确定后缀规则：若器官为 "jiechang"，后缀设置为 "jie"，其他器官后缀为器官名称
  suffix <- if(organ == "jiechang") "jie" else organ
  
  ## 处理治疗组样本：治疗组可能包含多行数据，合并后所有样本名称
  treat_samples_raw <- treat_data[[organ]]
  # 去除NA以及空字符串
  treat_samples_raw <- treat_samples_raw[!is.na(treat_samples_raw) & treat_samples_raw != ""]
  # 拆分每个单元格中以逗号分隔的多个样本
  treat_samples_list <- strsplit(as.character(treat_samples_raw), ", ")
  treat_samples <- unlist(treat_samples_list)
  # 对非血液样本添加后缀（血液样本特殊处理）
  if (organ != "blood") {
    treat_samples <- paste0(treat_samples, ".", suffix)
  } else {
    # 如果器官是blood, 假设只有治疗组需要添加后缀
    treat_samples <- paste0(treat_samples, ".", suffix)
  }
  
  ## 处理对照组样本
  control_samples_raw <- control_data[[organ]]
  control_samples_raw <- control_samples_raw[!is.na(control_samples_raw) & control_samples_raw != ""]
  control_samples_list <- strsplit(as.character(control_samples_raw), ", ")
  control_samples <- unlist(control_samples_list)
  # 对于blood以外的器官，对照组也添加后缀；若器官为blood，则直接使用对照组原样
  if (organ != "blood") {
    control_samples <- paste0(control_samples, ".", suffix)
  }
  
  # 检查样本是否都在 merged_data 列中
  all_samples <- c(treat_samples, control_samples)
  missing_samples <- all_samples[!all_samples %in% colnames(merged_data)]
  if(length(missing_samples) > 0) {
    warning(paste0("在器官 ", organ, " 中以下样本未在数据中找到: ", paste(missing_samples, collapse = ", ")))
    next
  }
  
  # 从 merged_data 中提取对应的表达矩阵
  expr_matrix <- merged_data[, all_samples]
  
  # 计算处理组和对照组的均值
  treat_means   <- rowMeans(expr_matrix[, treat_samples, drop = FALSE])
  control_means <- rowMeans(expr_matrix[, control_samples, drop = FALSE])
  
  # 筛选两组均值均大于 0 的基因（以免 log2(0) 出问题）
  keep_genes <- treat_means > 0 & control_means > 0
  if (sum(keep_genes) == 0) {
    warning(paste0("器官 ", organ, " 中没有符合条件的基因"))
    next
  }
  expr_matrix <- expr_matrix[keep_genes, ]
  gene_list   <- merged_data$SYMBOL[keep_genes]
  treat_means   <- treat_means[keep_genes]
  control_means <- control_means[keep_genes]
  
  # 计算log2FC：log2(处理组均值) - log2(对照组均值)
  log2fc <- log2(treat_means) - log2(control_means)
  
  # 整理结果数据框
  res_df <- data.frame(
    SYMBOL = gene_list,
    treat_mean = treat_means,
    control_mean = control_means,
    log2FC = log2fc,
    stringsAsFactors = FALSE
  )
  
  # 修改保存文件名的规则：
  # 如果器官为 "blood"，则文件名设为 "all_drugs_vs_PBS_limma_Blood_results.csv"
  # 否则，文件名设为 "all_drugs_vs_PBS_limma_[organ]_results.csv"
  output_file <- if(organ == "blood") {
    "all_drugs_vs_PBS_limma_Blood_results.csv"
  } else {
    paste0("all_drugs_vs_PBS_limma_", organ, "_results.csv")
  }
  
  # 保存结果到 CSV 文件
  write.csv(res_df, file = output_file, quote = FALSE, row.names = FALSE)
  
  cat("完成器官", organ, "的分析，结果保存在文件：", output_file, "\n")
}