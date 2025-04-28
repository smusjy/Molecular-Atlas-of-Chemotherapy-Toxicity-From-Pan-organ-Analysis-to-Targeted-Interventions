# 载入必要的包
library(readr)
library(dplyr)
library(readxl)

# 读取生化指标数据
Data_Biochemical <- read_excel("Data_Blood.xlsx")
cat("已读取生化指标数据\n")

# 设定ssGSEA结果的主文件夹路径
main_path <- "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results"

# 获取所有器官文件夹
organs <- list.dirs(main_path, full.names = FALSE, recursive = FALSE)
cat("发现", length(organs), "个器官文件夹\n")

# 创建一个列表来存储所有药物的数据
all_drug_data <- list()

# 获取所有药物名称
drug_organ <- read_csv("drug_organ.csv")
all_drugs <- unique(drug_organ$Drug)
cat("发现", length(all_drugs), "种药物\n")

# 提取样本编号的函数
extract_number <- function(x) {
  as.numeric(gsub("[^0-9]", "", x))
}

# 对每个药物进行处理
cat("\n开始处理每个药物的数据：\n")
for(drug_idx in seq_along(all_drugs)) {
  drug <- all_drugs[drug_idx]
  cat(sprintf("\n处理药物 %d/%d: %s\n", drug_idx, length(all_drugs), drug))
  
  # 创建列表存储当前药物在各器官的ssGSEA分数
  drug_organ_scores <- list()
  
  # 对每个器官进行处理
  cat("读取各器官数据:\n")
  for(organ_idx in seq_along(organs)) {
    organ <- organs[organ_idx]
    cat(sprintf("  %d/%d: %s ", organ_idx, length(organs), organ))
    
    if(drug == "PBS") {
      file_path <- file.path(main_path, organ, "PBS_merged_gsvascore.csv")
    } else {
      file_path <- file.path(main_path, organ, paste0(drug, "_merged_gsvascore.csv"))
    }
    
    # 检查文件是否存在
    if(file.exists(file_path)) {
      # 读取文件
      organ_data <- read.csv(file_path, row.names = 1)
      drug_organ_scores[[organ]] <- organ_data
      cat("✓\n")
    } else {
      cat("×\n")
    }
  }
  
  # 如果没有找到任何数据，跳过这个药物
  if(length(drug_organ_scores) == 0) {
    cat("未找到数据，跳过\n")
    next
  }
  
  cat("计算平均分数...\n")
  
  # 获取所有通路名称
  pathways <- rownames(drug_organ_scores[[1]])
  
  # 获取所有样本编号并排序
  all_samples <- unique(unlist(lapply(drug_organ_scores, colnames)))
  all_sample_numbers <- sapply(all_samples, extract_number)
  sorted_samples <- all_samples[order(all_sample_numbers)]
  
  # 创建一个矩阵来存储平均分数
  avg_scores <- matrix(NA, nrow = length(pathways), ncol = length(sorted_samples))
  rownames(avg_scores) <- pathways
  colnames(avg_scores) <- sorted_samples
  
  # 计算每个通路在所有器官的平均分数
  total_pathways <- length(pathways)
  for(p_idx in seq_along(pathways)) {
    if(p_idx %% 100 == 0) {
      cat(sprintf("  处理通路进度: %d/%d\n", p_idx, total_pathways))
    }
    
    pathway <- pathways[p_idx]
    for(sample in sorted_samples) {
      # 收集该样本在各个器官的得分
      sample_scores <- numeric()
      for(organ in names(drug_organ_scores)) {
        if(sample %in% colnames(drug_organ_scores[[organ]])) {
          sample_scores <- c(sample_scores, 
                             drug_organ_scores[[organ]][pathway, sample])
        }
      }
      # 如果有分数，计算平均值
      if(length(sample_scores) > 0) {
        avg_scores[pathway, sample] <- mean(sample_scores, na.rm = TRUE)
      }
    }
  }
  
  cat("处理样本匹配...\n")
  
  # 移除全NA的列
  avg_scores <- avg_scores[, colSums(is.na(avg_scores)) < nrow(avg_scores)]
  
  # 转换为数据框并添加通路名称列
  avg_scores_df <- as.data.frame(avg_scores)
  avg_scores_df$Pathway <- rownames(avg_scores_df)
  
  # 处理生化指标数据
  if(drug == "PBS") {
    drug_bio <- Data_Biochemical[Data_Biochemical$Drug == "PBS",]
    # 处理PBS样本名称对应关系
    drug_samples_bio <- drug_bio$Sample
    drug_samples_data <- colnames(avg_scores)
  } else {
    drug_bio <- Data_Biochemical[Data_Biochemical$Drug == drug,]
    # 处理药物样本对应关系
    drug_samples_bio <- drug_bio$Sample
    drug_samples_data <- colnames(avg_scores)
  }
  
  # 按数字大小排序
  drug_samples_bio <- drug_samples_bio[order(sapply(drug_samples_bio, extract_number))]
  drug_samples_data <- drug_samples_data[order(sapply(drug_samples_data, extract_number))]
  
  # 更新对应关系
  for(i in 1:length(drug_samples_bio)) {
    drug_bio$Sample[drug_bio$Sample == drug_samples_bio[i]] <- drug_samples_data[i]
  }
  
  # 存储结果
  all_drug_data[[drug]] <- list(
    scores = avg_scores_df,
    biochem = drug_bio
  )
}

# 计算相关性
cat("\n\n开始计算相关性...\n")
all_correlations <- list()

total_drugs <- length(names(all_drug_data))
for(drug_idx in seq_along(names(all_drug_data))) {
  drug <- names(all_drug_data)[drug_idx]
  if(drug == "PBS") next  # 跳过单独的PBS处理
  
  cat(sprintf("\n处理药物 %d/%d: %s\n", drug_idx, total_drugs, drug))
  
  drug_data <- all_drug_data[[drug]]
  pbs_data <- all_drug_data[["PBS"]]
  
  # 合并当前药物和PBS的数据
  combined_scores <- drug_data$scores
  combined_biochem <- rbind(drug_data$biochem, pbs_data$biochem)
  
  # 获取通路和生化指标
  pathways <- combined_scores$Pathway
  biochem_indices <- colnames(combined_biochem)[3:ncol(combined_biochem)]
  
  total_pathways <- length(pathways)
  for(p_idx in seq_along(pathways)) {
    if(p_idx %% 100 == 0) {
      cat(sprintf("  处理通路进度: %d/%d\n", p_idx, total_pathways))
    }
    
    pathway <- pathways[p_idx]
    pathway_scores <- as.numeric(combined_scores[combined_scores$Pathway == pathway, 
                                                 1:(ncol(combined_scores)-1)])
    
    for(bio in biochem_indices) {
      bio_values <- as.numeric(combined_biochem[[bio]])
      
      # 移除NA值
      valid_indices <- !is.na(pathway_scores) & !is.na(bio_values)
      
      # 尝试计算相关性，如果失败则填充NA
      tryCatch({
        if(sum(valid_indices) >= 3) {  # 确保至少有3个有效配对
          cor_test <- cor.test(pathway_scores[valid_indices], 
                               bio_values[valid_indices], 
                               method = "spearman", 
                               exact = FALSE)
          
          all_correlations[[length(all_correlations) + 1]] <- data.frame(
            Drug = drug,
            Pathway = pathway,
            Biochemical = bio,
            Correlation = round(cor_test$estimate, 3),
            P_value = round(cor_test$p.value, 3),
            Valid_pairs = sum(valid_indices)
          )
        } else {
          all_correlations[[length(all_correlations) + 1]] <- data.frame(
            Drug = drug,
            Pathway = pathway,
            Biochemical = bio,
            Correlation = NA,
            P_value = NA,
            Valid_pairs = sum(valid_indices)
          )
        }
      }, error = function(e) {
        all_correlations[[length(all_correlations) + 1]] <- data.frame(
          Drug = drug,
          Pathway = pathway,
          Biochemical = bio,
          Correlation = NA,
          P_value = NA,
          Valid_pairs = sum(valid_indices)
        )
      })
    }
  }
}

# 合并所有结果
cat("\n合并所有结果...\n")
final_results <- do.call(rbind, all_correlations)

# 保存结果
cat("保存结果...\n")
write.csv(final_results, "Drug_pathway_biochemical_correlations.csv", row.names = FALSE)
cat("\n分析完成！结果已保存到 Drug_pathway_biochemical_correlations.csv\n")