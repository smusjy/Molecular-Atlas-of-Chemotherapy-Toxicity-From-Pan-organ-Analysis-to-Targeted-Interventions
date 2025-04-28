# -------------------------------
# 1. 加载所需的包
# -------------------------------
library(limma)
library(readr)
library(tidyverse)

# -------------------------------
# 2. 设置全局参数
# -------------------------------
delta <- 0.1  # 安全边际，确保平移后的数据均 > 0

# 定义各器官信息（请根据实际情况调整）
# name: 分析结果中使用的器官名称
# dir: 数据存放的文件夹名称，例如 “肺”、“血”、“心” 等
# suffix: 用于文件名匹配（假定文件名中包含该字符串，如 "lung"）
organs <- list(
  list(name = "Lung",           dir = "肺",      suffix = "lung"),
  list(name = "Blood",          dir = "血",      suffix = "blood"),
  list(name = "Heart",          dir = "心",      suffix = "xin"),
  list(name = "BloodVessel",    dir = "血管",    suffix = "guan"),
  list(name = "Liver",          dir = "肝",      suffix = "gan"),
  list(name = "Testis",         dir = "睾丸",    suffix = "gao"),
  list(name = "Brain",          dir = "脑",      suffix = "nao"),
  list(name = "Bone",           dir = "骨",      suffix = "gu"),
  list(name = "Colon",          dir = "结肠",    suffix = "jiechang"),
  list(name = "SkeletalMuscle", dir = "肌肉",    suffix = "ji"),
  list(name = "Stomach",        dir = "胃",      suffix = "wei"),
  list(name = "Skin",           dir = "皮肤",    suffix = "pi"),
  list(name = "Prostate",       dir = "前列腺",  suffix = "qian"),
  list(name = "Kidney",         dir = "肾",      suffix = "shen")
)

# 假定所有器官数据均存放在 base_path 下各器官的文件夹中
base_path <- "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results"

# -------------------------------
# 3. 第一轮循环：读取每个器官数据，
#    通过通路名匹配后计算原始（未加平移）处理组和对照组分别的均值，
#    并将这些值存储到一个 data.frame（以后根据通路名来合并）
# -------------------------------
organ_data_list <- list()  # 存储每个器官的数据信息及对应的原始均值
all_scores <- c()         # 用于累计所有器官的数值计算全局最小值

for (org in organs) {
  
  cat("Processing organ:", org$name, "\n")
  
  # 当前器官的工作目录
  base_dir <- file.path(base_path, org$dir)
  
  ## 读取对照组数据：假设文件名为 "PBS_merged_gsvascore.csv"
  control_file <- list.files(path = base_dir, 
                             pattern = "^PBS_merged_gsvascore\\.csv$", 
                             full.names = TRUE)
  if (length(control_file) == 0) {
    stop(paste("Control file not found in", org$name, "directory!"))
  }
  control_data <- read.csv(control_file[1], check.names = FALSE)
  rownames(control_data) <- control_data[,1]
  control_data <- control_data[,-1]
  
  ## 读取处理组数据：忽略以 "PBS_" 开头的文件，其余文件格式为 "DrugName_merged_gsvascore.csv"
  all_files <- list.files(path = base_dir, 
                          pattern = "merged_gsvascore\\.csv$", 
                          full.names = TRUE)
  
  treatment_files <- all_files[ !grepl("^PBS_", basename(all_files)) ]
  if (length(treatment_files) == 0) {
    stop(paste("No treatment files found in", org$name, "directory!"))
  }
  
  # 获取 drug 名称，假定 drug 名称为文件名第一个 "_" 之前的部分
  drug_names <- sapply(basename(treatment_files), function(x){
    strsplit(x, "_")[[1]][1]
  })
  unique_drugs <- unique(drug_names)
  
  # 将每个 drug 对应的文件数据按列合并
  all_data <- NULL   # 存储所有药物（处理组）的数据
  for(drug in unique_drugs) {
    pattern_drug <- paste0("^", drug, "_merged_gsvascore\\.csv$")
    drug_file <- list.files(path = base_dir, pattern = pattern_drug, full.names = TRUE)
    if(length(drug_file) > 0) {
      drug_data <- read.csv(drug_file[1], check.names = FALSE)
      rownames(drug_data) <- drug_data[, 1]
      drug_data <- drug_data[, -1]
      if(is.null(all_data)){
        all_data <- drug_data
      } else {
        # 按照通路名称匹配，再合并
        common_paths <- intersect(rownames(all_data), rownames(drug_data))
        if(length(common_paths)==0){
          stop(paste("No common pathways found for drug", drug, "in", org$name))
        }
        # 保持相同顺序：排序后利用 match
        common_paths <- sort(common_paths)
        all_data <- all_data[match(common_paths, rownames(all_data)), , drop = FALSE]
        drug_data <- drug_data[match(common_paths, rownames(drug_data)), , drop = FALSE]
        all_data <- cbind(all_data, drug_data)
      }
    }
  }
  
  # 利用通路名称匹配处理组与对照组：提取共同的通路，并保证顺序一致
  common_paths <- intersect(rownames(all_data), rownames(control_data))
  if(length(common_paths)==0){
    stop(paste("No common pathways found between treatment and control in", org$name))
  }
  common_paths <- sort(common_paths)
  all_data <- all_data[match(common_paths, rownames(all_data)), , drop = FALSE]
  control_data <- control_data[match(common_paths, rownames(control_data)), , drop = FALSE]
  
  # 合并前记录原始均值
  original_treatment_mean <- rowMeans(all_data)
  original_control_mean <- rowMeans(control_data)
  
  # 保存“原始均值”结果为一个 data.frame（以通路名称为关键字）
  original_means_df <- data.frame(
    Pathway = common_paths,
    original_treatment_mean = original_treatment_mean,
    original_control_mean = original_control_mean,
    stringsAsFactors = FALSE
  )
  
  # 合并处理组数据和对照组数据，保证行是共同的通路
  final_data <- cbind(all_data, control_data)
  
  # 将当前器官的所有数据保存至列表中
  organ_data_list[[org$name]] <- list(
    final_data = final_data,                 # 未加 k 调整的数据
    all_data = all_data,                     # 处理组数据
    control_data = control_data,             # 对照组数据
    original_means_df = original_means_df,     # 以通路为 key 储存的均值
    base_dir = base_dir
  )
  
  # 累计所有数值（用于后续全局确定 k）
  all_scores <- c(all_scores, as.numeric(as.matrix(final_data)))
}

# 计算全局最小值并确定平移常数 k
min_val <- min(all_scores, na.rm = TRUE)
if(min_val <= 0){
  k <- -min_val + delta
} else {
  k <- 0
}
cat("Global min value =", min_val, "\n")
cat("Determined constant k =", k, "\n")

# -------------------------------
# 4. 第二轮循环：对每个器官进行 limma 分析
#    对加 k 后的数据计算 log₂FC，然后根据通路名称将之前保存的原始均值合并进结果
# -------------------------------
for(org_name in names(organ_data_list)){
  
  cat("Performing limma analysis for organ:", org_name, "\n")
  
  # 获取当前器官数据
  data_list <- organ_data_list[[org_name]]
  final_data <- data_list$final_data   # 原始未调整数据（按通路名称排序一致）
  all_data <- data_list$all_data       # 处理组数据（匹配后的）
  control_data <- data_list$control_data  # 对照组数据
  original_means_df <- data_list$original_means_df  # 含每个通路原始均值的 data.frame
  
  # 调整数据：加上常数 k，使得数据均 > 0（便于对数转换）
  adjusted_final_data <- final_data + k
  
  # 计算调整后数据中处理组与对照组的均值（此处依然基于通路名称匹配，行顺序统一）
  treatment_mean_adjusted <- rowMeans(adjusted_final_data[, colnames(all_data), drop = FALSE])
  control_mean_adjusted <- rowMeans(adjusted_final_data[, colnames(control_data), drop = FALSE])
  
  # 计算 log₂FC（基于调整后的数据均值）
  log2FC_adjusted <- log2(treatment_mean_adjusted / control_mean_adjusted)
  
  # 构造分组向量：前半部分是处理组，后半部分是对照组
  group <- factor(c(rep("treatment", ncol(all_data)),
                    rep("control", ncol(control_data))))
  
  # 建立设计矩阵并使用 limma 分析
  design <- model.matrix(~ group)
  fit <- lmFit(adjusted_final_data, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  
  # 为了便于 merge，将结果中的行名（通路名）转换为一列
  results_df <- tibble::rownames_to_column(results, var="Pathway")
  
  # 生成手动计算的 log₂FC 结果，并生成一个 data.frame
  manual_log2FC_df <- data.frame(Pathway = rownames(final_data),
                                 manual_log2FC = log2FC_adjusted,
                                 stringsAsFactors = FALSE)
  
  # 基于通路名称将原始均值与手动 log₂FC 进行合并
  merged_df <- merge(manual_log2FC_df, original_means_df, by = "Pathway", all.x = TRUE)
  
  # 将合并后的信息再与 limma 分析结果合并
  final_results <- merge(results_df, merged_df, by = "Pathway", all.x = TRUE)
  
  # 如果需要，可以调整 final_results 排序或行名
  rownames(final_results) <- final_results$Pathway
  
  # 输出结果——请根据实际需要调整保存路径
  output_file <- file.path("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/Figure 2F",
                           paste0("all_drugs_vs_PBS_limma_", org_name, "_results.csv"))
  write.csv(final_results, output_file, row.names = TRUE)
  cat("Results saved to:", output_file, "\n")
}

cat("All organ analyses completed.\n")