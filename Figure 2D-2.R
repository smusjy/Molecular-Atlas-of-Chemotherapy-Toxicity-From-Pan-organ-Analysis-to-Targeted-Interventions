library(dplyr)
library(tidyr)
library(ggplot2)

# 设置工作目录
setwd("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/Figure 2D/log2FC")

# 定义器官拼音与英文名称的映射字典
organ_true_names <- c(
  "xin"   = "Heart",
  "gan"   = "Liver",
  "fei"   = "Lung",
  "shen"  = "Kidney",
  "pi"    = "Skin",
  "gu"    = "Bone",
  "guan"  = "BloodVessel",
  "ji"    = "SkeletalMuscle",
  "jiechang"    = "Colon",
  "nao"    = "Brain",
  "qian"    = "Prostate",
  "wei"    = "Stomach",
  "gao" = "Testis",
  "Blood" = "Blood"
)

# 获取所有csv文件
files <- list.files(pattern = ".*limma.*\\.csv$")

# 创建列表存储所有原始数据
all_data_list <- list()

# 读取所有数据，并将器官名字替换为英文
for(file in files) {
  # 从文件名中提取器官的拼音名称
  organ <- gsub(".*limma_(.+)_results\\.csv", "\\1", file)
  # 如果映射字典中存在，则替换为英文名称
  if(organ %in% names(organ_true_names)) {
    organ <- organ_true_names[organ]
  }
  data <- read.csv(file, header = TRUE)
  colnames(data)[1] <- "gene"
  # 过滤掉不需要的基因
  data <- data %>%
    filter(!grepl("Rik$", gene),  # 过滤掉Rik结尾的基因
           !grepl("^Gm", gene))   # 过滤掉Gm开头的基因
  all_data_list[[organ]] <- data
}

# 修改函数以考虑p值显著性
calculate_zscores_and_pvalues <- function(target_organ, all_data_list, p_threshold = 0.05) {
  # 获取目标器官的数据
  target_data <- all_data_list[[target_organ]]
  
  # 创建结果数据框
  results <- data.frame(
    gene = target_data$gene,
    logFC = target_data$logFC
  )
  
  # 计算每个基因的z-score和p值
  for(gene in results$gene) {
    all_logFCs <- numeric()
    target_logFC <- target_data$logFC[target_data$gene == gene]
    
    # 收集所有器官的logFC值
    for(organ in names(all_data_list)) {
      organ_data <- all_data_list[[organ]]
      if(gene %in% organ_data$gene) {
        all_logFCs <- c(all_logFCs, 
                        organ_data$logFC[organ_data$gene == gene])
      }
    }
    
    # 计算z-score
    if(length(all_logFCs) > 1) {
      mean_logFC <- mean(all_logFCs, na.rm = TRUE)
      sd_logFC <- sd(all_logFCs, na.rm = TRUE)
      
      if(!is.na(sd_logFC) && sd_logFC != 0) {
        results$zscore[results$gene == gene] <- 
          (target_logFC - mean_logFC) / sd_logFC
      } else {
        results$zscore[results$gene == gene] <- NA
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
          results$pvalue[results$gene == gene] <- test_result$p.value
        }, error = function(e) {
          results$pvalue[results$gene == gene] <- NA
        })
      } else {
        results$pvalue[results$gene == gene] <- NA
      }
    } else {
      results$zscore[results$gene == gene] <- NA
      results$pvalue[results$gene == gene] <- NA
    }
  }
  
  # 获取top 5基因（仅考虑z-score，过滤掉不显著的）
  top_genes <- results %>%
    filter(!is.na(zscore), !is.na(pvalue), pvalue < p_threshold) %>%
    arrange(desc(zscore)) %>% 
    head(5) %>%
    pull(gene)
  
  return(list(results = results, top_genes = top_genes))
}

# 对每个器官计算z-score和p值，并获取top基因
all_top_genes <- c()
organ_results <- list()

for(organ in names(all_data_list)) {
  result <- calculate_zscores_and_pvalues(organ, all_data_list, p_threshold = 0.05)
  organ_results[[organ]] <- result$results
  all_top_genes <- c(all_top_genes, result$top_genes)
}

# 为重复基因添加编号
gene_counts <- table(all_top_genes)
genes_with_suffix <- all_top_genes
counter <- list()

for(gene in names(gene_counts)) {
  if(gene_counts[gene] > 1) {
    counter[[gene]] <- 1
    indices <- which(all_top_genes == gene)
    for(idx in indices) {
      genes_with_suffix[idx] <- paste0(gene, "_", counter[[gene]])
      counter[[gene]] <- counter[[gene]] + 1
    }
  }
}

all_top_genes <- genes_with_suffix

# 创建最终用于可视化的数据框
final_data <- data.frame()

for(organ in names(all_data_list)) {
  organ_data <- organ_results[[organ]]
  
  for(i in seq_along(all_top_genes)) {
    gene <- all_top_genes[i]
    original_gene <- sub("_\\d+$", "", gene)
    
    if(original_gene %in% organ_data$gene) {
      final_data <- rbind(final_data, data.frame(
        organ = organ,
        gene = gene,
        zscore = organ_data$zscore[organ_data$gene == original_gene],
        pvalue = organ_data$pvalue[organ_data$gene == original_gene]
      ))
    }
  }
}

# 按照器官英文名称首字母排序
unique_organs <- unique(final_data$organ)
organ_order <- unique_organs[order(sapply(unique_organs, function(x) { substr(x, 1, 1) }))]

# 更新最终数据中的器官因子水平
final_data$organ <- factor(final_data$organ, levels = organ_order)

# 新的基因排序方法
gene_order <- c()
for(organ in rev(organ_order)) {
  organ_genes <- final_data %>%
    filter(organ == !!organ) %>%
    group_by(gene) %>%
    slice_max(order_by = zscore, n = 1) %>%
    arrange(desc(zscore)) %>%
    pull(gene) %>%
    head(5)
  
  gene_order <- c(organ_genes, gene_order)
}

# 更新final_data中gene的因子水平
final_data$gene <- factor(final_data$gene, levels = gene_order)

# 创建完整的可视化矩阵
# Z-score矩阵
zscore_matrix <- final_data %>%
  select(organ, gene, zscore) %>%
  pivot_wider(
    names_from = gene,
    values_from = zscore
  )

# p-value矩阵
pvalue_matrix <- final_data %>%
  select(organ, gene, pvalue) %>%
  pivot_wider(
    names_from = gene,
    values_from = pvalue
  )

# 保存矩阵
write.csv(zscore_matrix, "zscore_matrix.csv", row.names = FALSE)
write.csv(pvalue_matrix, "pvalue_matrix.csv", row.names = FALSE)

# 创建带标签的可视化图
ggplot(final_data, aes(x = gene, y = factor(organ, levels = rev(levels(final_data$organ))))) +
  geom_point(aes(color = zscore, 
                 size = -log10(pvalue)),
             na.rm = TRUE) +
  scale_color_gradient2(
    low = "#204A8D",
    mid = "white",
    high = "#e31a1c",
    midpoint = 0,
    name = "Z-score"
  ) +
  scale_size_continuous(
    name = "-log10(p-value)",
    range = c(0.5, 5)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Gene", 
       y = "Organ",
       title = "Gene Expression Z-scores Across Organs")

# 保存带标签的图片
ggsave("organ_specific_expression_zscores.pdf", width = 14, height = 8, dpi = 300)

# 创建无标签版本的可视化图
ggplot(final_data, aes(x = gene, y = factor(organ, levels = rev(levels(final_data$organ))))) +
  geom_point(aes(color = zscore, 
                 size = -log10(pvalue)),
             na.rm = TRUE) +
  scale_color_gradient2(
    low = "#204A8D",
    mid = "white",
    high = "#e31a1c",
    midpoint = 0,
    name = "Z-score"
  ) +
  scale_size_continuous(
    name = "-log10(p-value)",
    range = c(0.5, 5)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "Gene", 
       y = "",
       title = "Gene Expression Z-scores Across Organs")

# 保存无标签版本图片
ggsave("organ_specific_expression_zscores_no_labels.pdf", width = 12, height = 3.3, dpi = 300)

# 输出带有z分数的TOP5基因列表
organ_top5_with_zscore <- final_data %>%
  group_by(organ) %>%
  slice_max(order_by = zscore, n = 5) %>%
  select(organ, gene, zscore) %>%
  arrange(organ)

# 保存TOP5基因列表
write.csv(organ_top5_with_zscore, "organ_specific_top5_genes_with_zscore.csv", row.names = FALSE)