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


# PART 2: 泛器官水平比较化疗药物处理组与PBS对照组的 log2FC #################
library(tidyr)

# 读入药物-器官分组表（文件内容见上文，列名依次为：Drug, blood, xin, guan, fei, gan, gao, nao, gu, jiechang, ji, wei, pi, qian, shen）
drug_organ <- read_csv("drug_organ.csv")

# 获取除 Drug 外的所有器官列
organ_cols <- names(drug_organ)[-1]

# 分离对照组（PBS）和处理组（非PBS）
treat_data <- drug_organ %>% filter(Drug != "PBS")
control_data <- drug_organ %>% filter(Drug == "PBS")

## 构建泛器官水平下的对照组样本（PBS）
control_samples_all <- c()
for(organ in organ_cols) {
  control_samples_raw <- control_data[[organ]]
  # 去除NA和空字符串
  control_samples_raw <- control_samples_raw[!is.na(control_samples_raw) & control_samples_raw != ""]
  if(length(control_samples_raw) == 0) next
  
  # 拆分单元格中以逗号及空格分隔的多个样本
  control_samples_list <- strsplit(as.character(control_samples_raw), ", ")
  control_samples <- unlist(control_samples_list)
  
  # 对于所有器官添加后缀：  
  # 如果器官为blood，检查PBS样本是否已以 "EAP.PBS" 开头，如是则保留原样，否则添加".blood"
  if(organ == "blood") {
    control_samples <- sapply(control_samples, function(x) {
      if(grepl("^EAP\\.PBS", x)) { 
        x 
      } else { 
        paste0(x, ".blood") 
      }
    })
    control_samples <- as.character(control_samples)
  } else {
    # 如果器官为jiechang，用 "jie"，否则用器官名本身作后缀
    suffix <- if(organ == "jiechang") "jie" else organ
    control_samples <- paste0(control_samples, ".", suffix)
  }
  control_samples_all <- c(control_samples_all, control_samples)
}
control_samples_all <- unique(control_samples_all)

## 针对每个化疗药物（非PBS）进行泛器官水平分析
for(i in 1:nrow(treat_data)) {
  drug <- treat_data$Drug[i]
  cat("开始分析药物:", drug, "\n")
  
  # 针对当前药物，遍历所有器官列，合并所有治疗组样本
  treat_samples_all <- c()
  for(organ in organ_cols) {
    treat_samples_raw <- treat_data[i, organ]  # 当前药物在该器官列的样本信息
    treat_samples_raw <- as.character(treat_samples_raw)
    if(is.na(treat_samples_raw) || treat_samples_raw == "") next
    
    # 拆分以逗号和空格分隔的多个样本
    treat_samples_list <- strsplit(treat_samples_raw, ", ")
    treat_samples <- unlist(treat_samples_list)
    
    # 对所有器官添加后缀：对于blood统一添加".blood"，对于jiechang添加"jie"，其他器官直接使用器官名作后缀
    if (organ == "blood") {
      treat_samples <- paste0(treat_samples, ".blood")
    } else {
      suffix <- if(organ == "jiechang") "jie" else organ
      treat_samples <- paste0(treat_samples, ".", suffix)
    }
    treat_samples_all <- c(treat_samples_all, treat_samples)
  }
  # 确保样本名称唯一
  treat_samples_all <- unique(treat_samples_all)
  
  # 合并处理组与对照组所有样本列表
  all_samples <- c(treat_samples_all, control_samples_all)
  
  # 检查样本是否都在 merged_data 列中
  missing_samples <- all_samples[!all_samples %in% colnames(merged_data)]
  if(length(missing_samples) > 0) {
    warning(paste0("药物 ", drug, " 中以下样本未在数据中找到: ", paste(missing_samples, collapse = ", ")))
    # 保留存在的样本
    all_samples <- intersect(all_samples, colnames(merged_data))
  }
  if(length(all_samples) == 0) {
    warning(paste0("药物 ", drug, " 中没有任何匹配的样本，跳过当前药物"))
    next
  }
  
  # 从 merged_data 中提取对应样本的表达矩阵
  expr_matrix <- merged_data[, all_samples]
  cat("expr_matrix维度: ", paste(dim(expr_matrix), collapse = " x "), "\n")
  
  # 计算治疗组和对照组的均值
  treat_means   <- rowMeans(expr_matrix[, treat_samples_all, drop = FALSE])
  control_means <- rowMeans(expr_matrix[, control_samples_all, drop = FALSE])
  
  # 筛选治疗组和对照组均值大于 0 的基因（避免 log2(0) 问题）
  keep_genes <- treat_means > 0 & control_means > 0
  if (sum(keep_genes) == 0) {
    warning(paste0("药物 ", drug, " 中没有符合条件的基因"))
    next
  }
  expr_matrix <- expr_matrix[keep_genes, ]
  gene_list   <- merged_data$SYMBOL[keep_genes]
  treat_means   <- treat_means[keep_genes]
  control_means <- control_means[keep_genes]
  
  # 计算 log2FC：log2(治疗组均值) - log2(对照组均值)
  log2fc <- log2(treat_means) - log2(control_means)
  
  # 整理结果数据框
  res_df <- data.frame(
    SYMBOL = gene_list,
    treat_mean = treat_means,
    control_mean = control_means,
    log2FC = log2fc,
    stringsAsFactors = FALSE
  )
  
  # 保存文件命名规则：all_drugs_vs_PBS_limma_<drug>_results.csv
  output_file <- paste0("all_drugs_vs_PBS_limma_", drug, "_results.csv")
  
  # 保存结果到 CSV 文件
  write.csv(res_df, file = output_file, quote = FALSE, row.names = FALSE)
  
  cat("完成药物", drug, "的分析，结果保存在文件：", output_file, "\n")
}