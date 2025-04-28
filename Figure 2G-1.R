library(readr)
library(tidyverse)

# -------------------------------
# 1. 基础设置
# -------------------------------

# 设置基础目录，该目录下包含各器官子文件夹
base_dir <- "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results"

# 导入药物信息（例如 CSV 表格中含有 Drug 列，不包含 PBS）
drug_organ <- read_csv("drug_organ.csv")
# 注意：此处 drugs 只取除 PBS 之外的药物名称
drugs <- drug_organ$Drug[c(1:13,15:17)]

# 定义所有器官名称与对应的文件夹后缀（需确保与文件名中的匹配一致）
organs <- c("肺", "血", "心", "血管", "肝", "睾丸", "脑", "骨", "结肠", "肌肉", "胃", "皮肤", "前列腺", "肾")
organ_files <- c("lung", "blood", "xin", "guan", "gan", "gao", "nao", "gu", "jiechang", "ji", "wei", "pi", "qian", "shen")

# 定义安全边际 delta，及稍后用于计算平移常数 k
delta <- 0.1  

# -------------------------------
# 2. 针对每个药物进行泛器官分析
# -------------------------------
for (current_drug in drugs) {
  
  cat("Analyzing drug:", current_drug, "\n")
  
  # 初始化全局数据（存储跨器官合并后的数据）
  all_drug_data <- NULL     # 处理组数据（行：通路）
  all_control_data <- NULL  # 对照组数据
  
  # 同时记录每个器官中（原始数据）的均值信息，以便后续对照
  means_list <- list()
  
  # 针对各个器官
  for (i in seq_along(organs)) {
    organ <- organs[i]
    organ_file <- organ_files[i]
    
    # 当前器官目录
    current_dir <- file.path(base_dir, organ)
    
    # 读取当前器官下处理组数据：命名例如 "Gemcitabine_merged_gsvascore.csv"
    treatment_pattern <- paste0("^", current_drug, "_merged_gsvascore\\.csv$")
    treatment_files <- list.files(current_dir, pattern = treatment_pattern, full.names = TRUE)
    
    # 读取当前器官下对照组数据：命名例如 "PBS_merged_gsvascore.csv"
    control_pattern <- "^PBS_merged_gsvascore\\.csv$"
    control_files <- list.files(current_dir, pattern = control_pattern, full.names = TRUE)
    
    if(length(treatment_files) == 0 || length(control_files) == 0){
      cat("Warning: Missing treatment or control file in", organ, "\n")
      next
    }
    
    # 读取处理组和对照组数据
    drug_data_org <- read.csv(treatment_files[1], check.names = FALSE)
    rownames(drug_data_org) <- drug_data_org[,1]
    drug_data_org <- drug_data_org[,-1]
    
    control_data_org <- read.csv(control_files[1], check.names = FALSE)
    rownames(control_data_org) <- control_data_org[,1]
    control_data_org <- control_data_org[,-1]
    
    # 按通路名称取处理组和对照组的公共部分，保证数据行对应（并按照字母排序）
    common_paths <- intersect(rownames(drug_data_org), rownames(control_data_org))
    if (length(common_paths) == 0) {
      cat("Warning: No common pathways in", organ, "\n")
      next
    }
    common_paths <- sort(common_paths)
    drug_data_org <- drug_data_org[match(common_paths, rownames(drug_data_org)), , drop = FALSE]
    control_data_org <- control_data_org[match(common_paths, rownames(control_data_org)), , drop = FALSE]
    
    # 记录当前器官下的原始均值信息
    means_list[[organ]] <- tibble(
      Pathway = common_paths,
      Mean_Treatment = rowMeans(drug_data_org),
      Mean_Control = rowMeans(control_data_org)
    )
    
    # 跨器官合并：如果当前器官数据存在，则和前面器官数据取公共通路后 cbind
    if (is.null(all_drug_data)) {
      all_drug_data <- drug_data_org
    } else {
      common_paths_global <- intersect(rownames(all_drug_data), rownames(drug_data_org))
      if (length(common_paths_global) == 0) {
        stop(paste("No common pathways across organs for drug", current_drug))
      }
      common_paths_global <- sort(common_paths_global)
      all_drug_data <- all_drug_data[match(common_paths_global, rownames(all_drug_data)), , drop = FALSE]
      drug_data_org <- drug_data_org[match(common_paths_global, rownames(drug_data_org)), , drop = FALSE]
      all_drug_data <- cbind(all_drug_data, drug_data_org)
    }
    
    if (is.null(all_control_data)) {
      all_control_data <- control_data_org
    } else {
      common_paths_global <- intersect(rownames(all_control_data), rownames(control_data_org))
      if (length(common_paths_global) == 0) {
        stop(paste("No common pathways across organs for control in drug", current_drug))
      }
      common_paths_global <- sort(common_paths_global)
      all_control_data <- all_control_data[match(common_paths_global, rownames(all_control_data)), , drop = FALSE]
      control_data_org <- control_data_org[match(common_paths_global, rownames(control_data_org)), , drop = FALSE]
      all_control_data <- cbind(all_control_data, control_data_org)
    }
    
  }  # end for each organ
  
  # 在跨器官层面，取处理组和对照组数据的公共通路
  overall_common_paths <- intersect(rownames(all_drug_data), rownames(all_control_data))
  if(length(overall_common_paths) == 0){
    stop(paste("No common pathways found for drug", current_drug, "across all organs"))
  }
  overall_common_paths <- sort(overall_common_paths)
  all_drug_data <- all_drug_data[match(overall_common_paths, rownames(all_drug_data)), , drop = FALSE]
  all_control_data <- all_control_data[match(overall_common_paths, rownames(all_control_data)), , drop = FALSE]
  
  # -------------------------------
  # 新增：计算全局最小值 & 确定平移常数 k（保证数据加上 k 后均 > 0）
  # -------------------------------
  all_values <- c(as.numeric(as.matrix(all_drug_data)),
                  as.numeric(as.matrix(all_control_data)))
  min_val <- min(all_values, na.rm = TRUE)
  if (min_val <= 0) {
    k <- -min_val + delta
  } else {
    k <- 0
  }
  cat("For drug", current_drug, "global min value =", min_val, "\n")
  cat("Determined constant k =", k, "\n")
  
  # 对数据加上平移常数 k
  adjusted_drug_data <- all_drug_data + k
  adjusted_control_data <- all_control_data + k
  
  # -------------------------------
  # 计算跨器官中每个通路的全局均值（基于调整后的数据）
  # -------------------------------
  drug_means <- rowMeans(adjusted_drug_data)
  control_means <- rowMeans(adjusted_control_data)
  
  # 计算 log₂FC：log₂(Mean_Treatment_adjusted) - log₂(Mean_Control_adjusted)
  log2FC <- log2(drug_means) - log2(control_means)
  
  # 构建原始均值数据框（以通路为 key），这里采用全局匹配的均值（未做 k 调整，仅作为对照参考）
  original_means_df <- tibble(
    Pathway = overall_common_paths,
    Mean_Treatment = rowMeans(all_drug_data),
    Mean_Control = rowMeans(all_control_data)
  )
  
  # 构建手动计算的 log₂FC 表（基于 k 调整后的数据）
  manual_log2FC_df <- tibble(
    Pathway = overall_common_paths,
    Log2FC = log2FC
  )
  
  # 最终将 log₂FC 与原始均值数据根据通路名称合并
  final_results <- left_join(manual_log2FC_df, original_means_df, by = "Pathway")
  
  # 保存结果到指定目录
  output_file <- file.path("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/Figure 2F",
                           paste0(current_drug, "_vs_PBS_manual_log2FC.csv"))
  write.csv(final_results, output_file, row.names = FALSE)
  cat("Results for", current_drug, "saved to:", output_file, "\n")
}

cat("All drug analyses completed.\n")