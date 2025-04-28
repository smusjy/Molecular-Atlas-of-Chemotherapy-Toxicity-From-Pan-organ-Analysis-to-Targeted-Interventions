# 加载必要的包
library(tidyverse)
library(limma)
library(ggplot2)

# 设置工作目录
setwd("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 3/Results")

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

# 选择需要的列
all_data <- all_data[,c(1:24,28:29)]

# 创建空的数据框来存储每个器官的limma结果
organ_limma_results <- data.frame()

# 对每个器官进行limma分析
unique_organs <- unique(all_data$Organ)
cell_types <- colnames(all_data)[!colnames(all_data) %in% c("Organ", "Drug", "ID", "X")]

for(organ in unique_organs) {
  # 提取当前器官的数据
  organ_data <- all_data[all_data$Organ == organ,]
  
  # 创建分组信息（泛化疗 vs PBS）
  group <- ifelse(organ_data$Drug == "PBS", "PBS", "Chemo")
  
  # 对每个细胞类型进行差异分析
  for(cell in cell_types) {
    # 准备数据矩阵
    expr_data <- as.numeric(organ_data[[cell]])
    expr_matrix <- t(matrix(expr_data, ncol=1))
    
    # 创建设计矩阵
    design <- model.matrix(~factor(group))
    
    # 进行limma差异分析
    fit <- lmFit(expr_matrix, design)
    fit <- eBayes(fit)
    
    # 提取结果
    result <- topTable(fit, coef=2, number=Inf)
    
    # 存储结果
    organ_limma_results <- rbind(organ_limma_results, data.frame(
      Organ = organ,
      Cell_Type = cell,
      logFC = result$logFC
    ))
  }
}

# 计算每个器官的z分数
final_results <- data.frame()

for(organ in unique_organs) {
  for(cell in cell_types) {
    # 获取当前器官和细胞类型的logFC
    current_logFC <- organ_limma_results$logFC[organ_limma_results$Organ == organ & 
                                                 organ_limma_results$Cell_Type == cell]
    
    # 获取其他器官的logFC
    other_logFC <- organ_limma_results$logFC[organ_limma_results$Organ != organ & 
                                               organ_limma_results$Cell_Type == cell]
    
    # 计算z分数
    z_score <- (current_logFC - mean(other_logFC)) / sd(other_logFC)
    
    # 进行单样本Wilcoxon检验
    p_value <- wilcox.test(current_logFC, mu = mean(other_logFC), 
                           alternative = "two.sided")$p.value
    
    # 存储结果
    final_results <- rbind(final_results, data.frame(
      Organ = organ,
      Cell_Type = cell,
      Z_score = z_score,
      P_value = p_value
    ))
  }
}

# 调整p值
final_results$adj.P_value <- p.adjust(final_results$P_value, method = "BH")

# 创建气泡图
ggplot(final_results, aes(y = Organ, x = Cell_Type)) +
  geom_point(aes(size = -log10(P_value), color = Z_score)) +
  scale_color_gradient2(low = "#204A8D", mid = "white", high = "#e31a1c", midpoint = 0) +
  scale_size_continuous(range = c(0.5, 5)) +
  scale_x_discrete(labels = function(x) gsub("_CIBERSORT$|_", " ", x)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    y = "Organ",
    x = "Immune Cell Type",
    color = "Z-score",
    size = "-log10(P-value)"
  )

# 保存结果
ggsave("organ_specific_immune_changes.pdf", width = 12, height = 8)
write.csv(final_results, "organ_specific_immune_changes.csv", row.names = FALSE)