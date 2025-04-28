# 加载必要的包
library(tidyverse)
library(limma)
library(ggplot2)

# 设置工作目录
setwd("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 3/Results/Results")

# 获取所有文件列表
files <- list.files(pattern = "^.*_.*_cibersort.*\\.csv$") 

# 创建空的数据框来存储结果
all_data <- data.frame()

# 读取并合并所有文件
for(file in files) {
  info <- strsplit(file, "_")[[1]]
  organ <- info[1]
  drug <- info[2]
  
  temp_data <- read.csv(file, header = TRUE)
  temp_data$Organ <- organ
  temp_data$Drug <- drug
  
  all_data <- rbind(all_data, temp_data)
}

# 选择需要的列（免疫细胞列）
all_data <- all_data[,c(1:24,28:29)]

# 获取唯一的药物名称（除PBS外）
unique_drugs <- unique(all_data$Drug[all_data$Drug != "PBS"])
cell_types <- colnames(all_data)[3:24]

# 对每个药物进行limma分析
for(drug in unique_drugs) {
  # 创建空的数据框来存储当前药物的limma结果
  drug_limma_results <- data.frame()
  
  # 提取当前药物组和PBS组的数据
  drug_data <- all_data[all_data$Drug == drug,]
  pbs_data <- all_data[all_data$Drug == "PBS",]
  
  # 合并数据
  combined_data <- rbind(drug_data, pbs_data)
  
  # 对每个细胞类型进行差异分析
  for(cell in cell_types) {
    # 准备数据矩阵
    expr_data <- as.numeric(combined_data[[cell]])
    expr_matrix <- t(matrix(expr_data, ncol=1))
    
    # 创建分组信息
    group <- factor(c(rep("drug", nrow(drug_data)), 
                      rep("PBS", nrow(pbs_data))))
    
    # 创建设计矩阵
    design <- model.matrix(~0+group)
    colnames(design) <- levels(group)
    
    # 进行limma差异分析
    fit <- lmFit(expr_matrix, design)
    contrast.matrix <- makeContrasts(drug-PBS, levels=design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    
    # 提取结果
    result <- topTable(fit2, number=Inf)
    
    # 存储结果
    drug_limma_results <- rbind(drug_limma_results, data.frame(
      Cell_Type = cell,
      logFC = result$logFC,
      P.Value = result$P.Value,
      adj.P.Val = result$adj.P.Val,
      t = result$t,
      B = result$B
    ))
  }
  
  # 为每个药物保存单独的结果文件，使用与原始文件相似的命名模式
  write.csv(drug_limma_results, 
            paste0("Pan_", drug, "_cibersort_limma.csv"), 
            row.names = FALSE)
}