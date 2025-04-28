library(dplyr)
library(tidyr)
library(ggplot2)
library(dendextend)
library(ggtree)
library(aplot)
library(ape)
library(reshape2)

# 设置工作目录
setwd("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 3/TCR/Mm_TCR_data/Chemo_Data_tsv")

#1. 数据准备和处理 ############################################################
# 读取数据
organ_data <- read.csv("Organ_Diversity_log2FC.csv")
drug_data <- read.csv("Drug_Overall_Diversity_log2FC.csv")

# 计算Z-scores和p值的函数
calculate_zscores <- function(data, group_col, measure_cols) {
  results <- data
  
  for(col in measure_cols) {
    values <- data[[col]]
    mean_val <- mean(values, na.rm = TRUE)
    sd_val <- sd(values, na.rm = TRUE)
    
    # 计算Z-score
    results[[paste0(col, "_zscore")]] <- (values - mean_val) / sd_val
    
    # 计算p值（与平均值的t检验）
    results[[paste0(col, "_pvalue")]] <- sapply(values, function(x) {
      other_vals <- values[values != x]
      if(length(other_vals) >= 2) {
        tryCatch({
          t.test(other_vals, mu = x)$p.value
        }, error = function(e) NA)
      } else {
        NA
      }
    })
  }
  
  return(results)
}

# 处理器官数据
organ_results <- calculate_zscores(
  organ_data,
  "Organ",
  c("Shannon_log2FC", "Simpson_log2FC", "Chao1_log2FC")
)

# 处理药物数据
drug_results <- calculate_zscores(
  drug_data,
  "Drug",
  c("Shannon_log2FC", "Simpson_log2FC", "Chao1_log2FC")
)

#2. 数据转换为长格式 #########################################################
# 器官数据转换
organ_long <- organ_results %>%
  select(Organ, ends_with("_zscore"), ends_with("_pvalue")) %>%
  pivot_longer(
    cols = -Organ,
    names_to = c("Index", "Measure"),
    names_pattern = "(.+)_(zscore|pvalue)$",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Measure,
    values_from = Value
  ) %>%
  mutate(Index = gsub("_log2FC", "", Index))

# 药物数据转换
drug_long <- drug_results %>%
  select(Drug, ends_with("_zscore"), ends_with("_pvalue")) %>%
  pivot_longer(
    cols = -Drug,
    names_to = c("Index", "Measure"),
    names_pattern = "(.+)_(zscore|pvalue)$",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Measure,
    values_from = Value
  ) %>%
  mutate(Index = gsub("_log2FC", "", Index))

#3. 聚类准备 ##############################################################
# 器官数据矩阵准备
organ_matrix <- reshape2::dcast(organ_long, Organ~Index, value.var = "zscore")
rownames(organ_matrix) <- organ_matrix$Organ
organ_matrix <- organ_matrix[,-1]

# 药物数据矩阵准备
drug_matrix <- reshape2::dcast(drug_long, Drug~Index, value.var = "zscore")
rownames(drug_matrix) <- drug_matrix$Drug
drug_matrix <- drug_matrix[,-1]

# 执行层次聚类
organ_clust <- hclust(dist(organ_matrix))
index_clust <- hclust(dist(t(organ_matrix)))

drug_clust <- hclust(dist(drug_matrix))
drug_index_clust <- hclust(dist(t(drug_matrix)))

# 转换为phylo对象
organ_phylo <- as.phylo(organ_clust)
drug_phylo <- as.phylo(drug_clust)

#4. 创建聚类树 #############################################################
# 器官聚类树
p_organ <- ggtree(organ_phylo, branch.length="none") + 
  xlim(NA, 7) +
  theme_tree2()

# 药物聚类树
p_drug <- ggtree(drug_phylo, branch.length="none") + 
  xlim(NA, 7) +
  theme_tree2()

#5. 重新排序数据 ##########################################################
# 根据聚类顺序重新设置因子水平
organ_long$Organ <- factor(organ_long$Organ, 
                           levels = rownames(organ_matrix)[organ_clust$order])
organ_long$Index <- factor(organ_long$Index, 
                           levels = colnames(organ_matrix)[index_clust$order])

drug_long$Drug <- factor(drug_long$Drug, 
                         levels = rownames(drug_matrix)[drug_clust$order])
drug_long$Index <- factor(drug_long$Index, 
                          levels = colnames(drug_matrix)[drug_index_clust$order])

#6. 创建主图 ##############################################################
# 器官多样性气泡图
p_main_organ <- ggplot(organ_long, aes(x = Index, y = Organ)) +
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
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Diversity Index", y = "Organ")

# 药物多样性气泡图
p_main_drug <- ggplot(drug_long, aes(x = Index, y = Drug)) +
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
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Diversity Index", y = "Drug")

#7. 合并图形并保存 ########################################################
# 合并器官图
final_organ_plot <- p_main_organ %>%
  insert_left(p_organ, width = 0.2)

# 合并药物图
final_drug_plot <- p_main_drug %>%
  insert_left(p_drug, width = 0.2)

# 保存图形
ggsave("organ_diversity_clustered.pdf", final_organ_plot, width = 6, height = 5, dpi = 300)
ggsave("drug_diversity_clustered.pdf", final_drug_plot, width = 6, height = 5, dpi = 300)

# 保存聚类顺序信息
write.csv(data.frame(
  Organ = rownames(organ_matrix)[organ_clust$order],
  Cluster_Order = 1:length(organ_clust$order)
), "organ_clustering_order.csv", row.names = FALSE)

write.csv(data.frame(
  Drug = rownames(drug_matrix)[drug_clust$order],
  Cluster_Order = 1:length(drug_clust$order)
), "drug_clustering_order.csv", row.names = FALSE)

# 保存矩阵
write.csv(organ_matrix, "organ_zscore_matrix.csv")
write.csv(drug_matrix, "drug_zscore_matrix.csv")