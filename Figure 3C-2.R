library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtree)
library(aplot)
library(ape)
library(reshape2)

# 设置工作目录
setwd("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 3/Results/Panorgan")

# 获取所有csv文件
files <- list.files(pattern = ".*limma.*\\.csv$")

# 创建列表存储所有原始数据
all_data_list <- list()

# 读取所有数据
for(file in files) {
  organ <- gsub(".*limma_(.+)_results\\.csv", "\\1", file)
  data <- read.csv(file, header = TRUE)
  colnames(data)[1] <- "cell"
  all_data_list[[organ]] <- data
}

# 修改函数以使用单样本t检验
calculate_zscores_and_pvalues <- function(target_organ, all_data_list, p_threshold = 0.05) {
  # 获取目标器官的数据
  target_data <- all_data_list[[target_organ]]
  
  # 创建结果数据框
  results <- data.frame(
    cell = target_data$cell,
    logFC = target_data$logFC
  )
  
  # 计算每个细胞类型的z-score和p值
  for(cell in results$cell) {
    all_logFCs <- numeric()
    target_logFC <- target_data$logFC[target_data$cell == cell]
    
    # 收集所有器官的logFC值
    for(organ in names(all_data_list)) {
      organ_data <- all_data_list[[organ]]
      if(cell %in% organ_data$cell) {
        all_logFCs <- c(all_logFCs, 
                        organ_data$logFC[organ_data$cell == cell])
      }
    }
    
    # 计算z-score
    if(length(all_logFCs) > 1) {
      mean_logFC <- mean(all_logFCs, na.rm = TRUE)
      sd_logFC <- sd(all_logFCs, na.rm = TRUE)
      
      if(!is.na(sd_logFC) && sd_logFC != 0) {
        results$zscore[results$cell == cell] <- 
          (target_logFC - mean_logFC) / sd_logFC
      } else {
        results$zscore[results$cell == cell] <- NA
      }
      
      # 获取其他器官的logFC值（不包括目标器官）
      other_logFCs <- all_logFCs[all_logFCs != target_logFC]
      
      # 进行单样本t检验
      if(length(other_logFCs) >= 2) {
        tryCatch({
          test_result <- t.test(
            other_logFCs,
            mu = target_logFC,
            alternative = "two.sided"
          )
          results$pvalue[results$cell == cell] <- test_result$p.value
        }, error = function(e) {
          results$pvalue[results$cell == cell] <- NA
        })
      } else {
        results$pvalue[results$cell == cell] <- NA
      }
    } else {
      results$zscore[results$cell == cell] <- NA
      results$pvalue[results$cell == cell] <- NA
    }
  }
  
  return(list(results = results))
}

# 对每个器官计算z-score和p值
organ_results <- list()

for(organ in names(all_data_list)) {
  result <- calculate_zscores_and_pvalues(organ, all_data_list)
  organ_results[[organ]] <- result$results
}

# 创建最终可视化数据框
final_data <- data.frame()

# 获取所有unique的cells
all_cells <- unique(unlist(lapply(all_data_list, function(x) x$cell)))

for(organ in names(all_data_list)) {
  organ_data <- organ_results[[organ]]
  
  for(cell in all_cells) {
    if(cell %in% organ_data$cell) {
      final_data <- rbind(final_data, data.frame(
        organ = organ,
        cell = cell,
        zscore = organ_data$zscore[organ_data$cell == cell],
        pvalue = organ_data$pvalue[organ_data$cell == cell]
      ))
    }
  }
}

# 数据清理和处理
final_data$cell <- as.character(final_data$cell)
final_data$cell <- gsub("_CIBERSORT", "", final_data$cell)
final_data$cell <- gsub("_", " ", final_data$cell)

final_data$organ <- as.character(final_data$organ)
final_data$organ <- gsub("_cibersort_limma.csv", "", final_data$organ)
final_data$organ <- gsub("Pan_", "", final_data$organ)

# 移除特定细胞类型
final_data <- final_data[!final_data$cell %in% c("Mast cells activated", "T cells gamma delta"),]

# Prepare clustering data
df <- final_data[,c("organ", "cell", "zscore")]
df_matrix <- reshape2::dcast(df, organ~cell, value.var = "zscore")
rownames(df_matrix) <- df_matrix$organ
df_matrix <- df_matrix[,-1]

# Perform hierarchical clustering
organ_clust <- hclust(dist(df_matrix))
cell_clust <- hclust(dist(t(df_matrix)))

# 将聚类树转换为phylo对象
organ_phylo <- as.phylo(organ_clust)
cell_phylo <- as.phylo(cell_clust)

# Create dendrograms
p_organ <- ggtree(organ_phylo, branch.length="none") + 
  xlim(NA, 7) +
  theme_tree2()

p_cell <- ggtree(cell_phylo, branch.length="none") + 
  xlim(NA, 7) +
  layout_dendrogram()

# Reorder data based on clustering
final_data$organ <- factor(final_data$organ, 
                           levels = rownames(df_matrix)[organ_clust$order])
final_data$cell <- factor(final_data$cell, 
                          levels = colnames(df_matrix)[cell_clust$order])

# Main plot
p_main <- ggplot(final_data, 
                 aes(x = cell, 
                     y = organ)) +
  geom_point(aes(color = zscore,
                 size = -log10(pvalue)),
             na.rm = TRUE) +
  scale_color_gradient2(
    low = "#204A8D",
    mid = "white",
    high = "#e31a1c",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    breaks = seq(-2, 2, 1),
    labels = seq(-2, 2, 1),
    name = "Z-score",
    oob = scales::squish
  ) +
  scale_size_continuous(
    name = "-log10(p-value)",
    range = c(0.5, 5)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45,
                               hjust = 1,
                               size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Cell Type", 
       y = "Organ")

# 组合最终图形
final_plot <- p_main %>%
  insert_left(p_organ, width = 0.2) %>%
  insert_top(p_cell, height = 0.2)

# 保存图形
ggsave("organ_specific_expression_clustered.pdf", 
       final_plot, 
       width = 10, 
       height = 4.5, 
       dpi = 300)