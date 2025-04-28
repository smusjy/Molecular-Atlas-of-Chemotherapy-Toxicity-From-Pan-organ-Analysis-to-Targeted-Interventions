# -------------------------------
# 1. 加载所需的包
# -------------------------------
library(limma)
library(readxl)
library(tidyverse)

# -------------------------------
# 2. 设置全局参数及器官列表
# -------------------------------
delta <- 0.1  # 安全边际，保证加常数后所有数据均 > 0

# 定义每个器官的信息
# name: 英文名称，用于结果命名
# dir: 数据存放文件夹（下属文件夹名称，与 ssGSEA.Rdata 文件所在路径组成完整路径）
# col: 在 drug_organ.xlsx 中对应的列名，用于提取各组样本名（注意：excel 文件必须包含列 "Drug" 以及对应器官的列）
organs <- list(
  list(name = "Lung",           dir = "肺",      col = "fei"),
  list(name = "Kidney",         dir = "肾",      col = "shen"),
  list(name = "Brain",          dir = "脑",      col = "jiaozhi"),
  list(name = "SkeletalMuscle", dir = "肌肉",    col = "ji"),
  list(name = "Testis",         dir = "睾丸",    col = "gao"),
  list(name = "Stomach",        dir = "胃",      col = "wei"),
  list(name = "Liver",          dir = "肝",      col = "gan"),
  list(name = "Bone",           dir = "骨",      col = "gu"),
  list(name = "Heart",          dir = "心",      col = "xin"),
  list(name = "BloodVessel",    dir = "血管",    col = "neipi"),
  list(name = "Colon",          dir = "结肠",    col = "jie"),
  list(name = "Prostate",       dir = "前列腺",  col = "qian"),
  list(name = "Skin",           dir = "皮肤",    col = "biao")
)

# 数据文件所在的基础文件夹（ssGSEA 的数据均存放在各器官子文件夹中，每个文件夹下包含 ssGSEA.Rdata）
base_path <- "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 4/ssGSEA"

# -------------------------------
# 3. 第一轮循环：读取各器官数据并计算原始均值
# -------------------------------
organ_data_list <- list()  # 用于保存每个器官的数据及各通路原始均值
all_scores <- c()         # 用于累计所有器官中所有数值，用于确定全局最小值，进而确定平移常数 k

for (org in organs) {
  
  cat("Processing organ:", org$name, "\n")
  
  ## 构建 ssGSEA 数据完整路径，并加载数据
  rdata_file <- file.path(base_path, org$dir, "ssGSEA.Rdata")
  if (!file.exists(rdata_file)) {
    stop(paste("数据文件不存在：", rdata_file))
  }
  load(rdata_file)  # 加载后会生成变量 ssgsea
  
  ## 读取药物-器官分组表
  drug_organ <- read_xlsx("drug_organ.xlsx")
  # 当前器官提取两列：“Drug”和对应器官列（例如 "fei"）
  current_drug_organ <- drug_organ[, c("Drug", org$col)]
  
  ## 按原逻辑，处理组为前 13 行和 15-17 行，对照组为第14行
  treat_drugs_idx <- c(1:13, 15:17)
  control_idx <- 14
  
  treat <- current_drug_organ$Drug[treat_drugs_idx]
  treat_organ <- current_drug_organ[treat_drugs_idx, ]
  control_organ <- current_drug_organ[control_idx, ]
  
  ## 获取所有处理组样本名称（注意，每个单元格内样本名以“, ”分隔）
  all_treat_samples <- c()
  for (i in seq_along(treat)) {
    treat_samples <- unlist(strsplit(treat_organ[[org$col]][i], ", "))
    all_treat_samples <- c(all_treat_samples, treat_samples)
  }
  control_samples <- unlist(strsplit(control_organ[[org$col]], ", "))
  
  ## 构造当前器官分析所需的样本名称及数据
  current_samples <- c(all_treat_samples, control_samples)
  current_data <- ssgsea[, current_samples, drop = FALSE]
  
  ## 分离处理组和对照组数据，用于计算“原始均值”
  treatment_data <- current_data[, all_treat_samples, drop = FALSE]
  control_data   <- current_data[, control_samples, drop = FALSE]
  
  original_treatment_mean <- rowMeans(treatment_data)
  original_control_mean   <- rowMeans(control_data)
  
  original_means_df <- data.frame(
    Pathway = rownames(current_data),
    original_treatment_mean = original_treatment_mean,
    original_control_mean = original_control_mean,
    stringsAsFactors = FALSE
  )
  
  ## 保存当前器官分析所需的数据列表
  organ_data_list[[org$name]] <- list(
    final_data = current_data,         # 原始 ssGSEA 数据（行为通路）
    treatment_samples = all_treat_samples,
    control_samples   = control_samples,
    original_means_df = original_means_df
  )
  
  ## 累计该器官中的所有数值，用于后续确定全局最小值
  all_scores <- c(all_scores, as.numeric(as.matrix(current_data)))
}

## 计算所有器官数据的全局最小值，并确定平移常数 k（确保所有数值加 k 后均为正）
min_val <- min(all_scores, na.rm = TRUE)
if (min_val <= 0) {
  k <- -min_val + delta
} else {
  k <- 0
}
cat("Global min value =", min_val, "\n")
cat("Determined constant k =", k, "\n")

# -------------------------------
# 4. 第二轮循环：对每个器官进行 limma 分析
#    采用加 k 后的数据计算 log₂FC，同时将手动计算 LFC 与原始均值结果合并
# -------------------------------
# 保存结果的文件夹（请根据实际情况调整）
output_base <- "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 4/ssGSEA/ssgsea_comparison_limma_pandrug"

for (org_name in names(organ_data_list)) {
  
  cat("Performing limma analysis for organ:", org_name, "\n")
  
  data_list <- organ_data_list[[org_name]]
  final_data <- data_list$final_data
  treatment_samples <- data_list$treatment_samples
  control_samples <- data_list$control_samples
  original_means_df <- data_list$original_means_df
  
  ## 对数据加上常数 k，确保数据均为正，便于后续对数转换
  adjusted_data <- final_data + k
  
  ## 手动计算调整后的处理组和对照组均值，并据此计算 log2FC
  treatment_mean_adjusted <- rowMeans(adjusted_data[, treatment_samples, drop = FALSE])
  control_mean_adjusted   <- rowMeans(adjusted_data[, control_samples, drop = FALSE])
  log2FC_adjusted <- log2(treatment_mean_adjusted) - log2(control_mean_adjusted)
  
  ## 构建分组向量：处理组在前，对照组在后
  group <- factor(c(rep("treatment", length(treatment_samples)),
                    rep("control", length(control_samples))))
  
  ## 使用 limma 进行差异分析
  design <- model.matrix(~ group)
  fit <- lmFit(adjusted_data, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  results_df <- tibble::rownames_to_column(results, var = "Pathway")
  
  ## 整理手动计算的 log2FC 结果，转为 data.frame
  manual_log2FC_df <- data.frame(
    Pathway = rownames(final_data),
    manual_log2FC = log2FC_adjusted,
    stringsAsFactors = FALSE
  )
  
  ## 合并手动 log2FC 与原始均值信息（均基于未平移的数据计算）
  merged_df <- merge(manual_log2FC_df, original_means_df, by = "Pathway", all.x = TRUE)
  
  ## 将 limma 分析结果与上步信息合并
  final_results <- merge(results_df, merged_df, by = "Pathway", all.x = TRUE)
  rownames(final_results) <- final_results$Pathway
  
  ## 构造输出文件路径并保存结果
  output_file <- file.path(output_base, paste0("all_drugs_vs_PBS_limma_", org_name, "_results.csv"))
  write.csv(final_results, output_file, row.names = TRUE)
  cat("Results saved to:", output_file, "\n")
}

cat("All organ analyses completed.\n")