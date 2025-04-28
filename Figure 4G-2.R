library(dplyr)
library(tidyr) 
library(ggplot2)
library(readxl)

# 设置工作目录
setwd("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 4/ssGSEA/ssgsea_comparison_limma_panorgan/log2FC")

# 获取所有csv文件
files <- list.files(pattern = ".*limma.*\\.csv$")

# 创建列表存储所有原始数据 
all_data_list <- list()

# 读取所有数据
for(file in files) {
  organ <- gsub(".*limma_(.+)_results\\.csv", "\\1", file)
  data <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
  colnames(data)[1] <- "pathway"
  all_data_list[[organ]] <- data
}

# 读取选定的通路列表
selected_pathways <- read_excel("cell_panorgan_pathway.xlsx")
pathways_to_show <- unique(selected_pathways$pathway)

# 计算通路的z-score和p值
calculate_zscores_and_pvalues <- function(target_organ, all_data_list) {
  target_data <- all_data_list[[target_organ]]
  
  results <- data.frame(
    pathway = target_data$pathway,
    log2FC = target_data$log2FC
  )
  
  for(pathway in results$pathway) {
    all_log2FCs <- numeric()
    target_log2FC <- target_data$log2FC[target_data$pathway == pathway]
    
    for(organ in names(all_data_list)) {
      organ_data <- all_data_list[[organ]]
      if(pathway %in% organ_data$pathway) {
        all_log2FCs <- c(all_log2FCs, 
                        organ_data$log2FC[organ_data$pathway == pathway])
      }
    }
    
    if(length(all_log2FCs) > 1) {
      mean_log2FC <- mean(all_log2FCs, na.rm = TRUE)
      sd_log2FC <- sd(all_log2FCs, na.rm = TRUE)
      
      if(!is.na(sd_log2FC) && sd_log2FC != 0) {
        results$zscore[results$pathway == pathway] <- 
          (target_log2FC - mean_log2FC) / sd_log2FC
      } else {
        results$zscore[results$pathway == pathway] <- NA
      }
      
      other_log2FCs <- all_log2FCs[all_log2FCs != target_log2FC]
      
      if(length(other_log2FCs) > 0) {
        test_result <- wilcox.test(
          x = target_log2FC,
          y = other_log2FCs,
          alternative = "two.sided"
        )
        results$pvalue[results$pathway == pathway] <- test_result$p.value
      } else {
        results$pvalue[results$pathway == pathway] <- NA
      }
    } else {
      results$zscore[results$pathway == pathway] <- NA
      results$pvalue[results$pathway == pathway] <- NA
    }
  }
  
  return(results)
}

# 对每个器官计算z-score和p值
organ_results <- list()
final_data <- data.frame()

for(organ in names(all_data_list)) {
  results <- calculate_zscores_and_pvalues(organ, all_data_list)
  organ_results[[organ]] <- results
  
  # 只保留选定的通路
  for(pathway in pathways_to_show) {
    if(pathway %in% results$pathway) {
      final_data <- rbind(final_data, data.frame(
        organ = organ,
        pathway = pathway,
        zscore = results$zscore[results$pathway == pathway],
        pvalue = results$pvalue[results$pathway == pathway]
      ))
    }
  }
}

# 读取通路名映射文件
pathway_mapping <- read_excel("cell_panorgan_pathway_2.xlsx")
names(pathway_mapping) <- c("new_name", "old_name")

# 清理通路名称的函数
clean_pathway_names <- function(pathway) {
  lowercase_words <- c("a", "an", "the", "and", "but", "or", "for", "nor", 
                       "on", "at", "to", "by", "of", "in", "with", "within")
  
  if(!grepl("^(GO_|REACTOME_|GOBP_)", pathway)) {
    return(pathway)
  }
  
  cleaned <- gsub("^(GO_|REACTOME_|GOBP_)", "", pathway)
  cleaned <- tolower(cleaned)
  cleaned <- gsub("_", " ", cleaned)
  
  cleaned <- gsub("(^|\\s)(\\w)", "\\1\\U\\2", cleaned, perl = TRUE)
  
  words <- strsplit(cleaned, " ")[[1]]
  for(i in 2:length(words)) {
    if(tolower(words[i]) %in% lowercase_words) {
      words[i] <- tolower(words[i])
    }
  }
  cleaned <- paste(words, collapse = " ")
  
  return(cleaned)
}

# 修改后的映射函数，添加编号以处理重复项
map_pathway_names <- function(pathway) {
  matching_row <- pathway_mapping[pathway_mapping$old_name == pathway, ]
  if(nrow(matching_row) > 0) {
    return(matching_row$new_name[1])
  } else {
    cleaned <- clean_pathway_names(pathway)
    # 为特定的重复通路添加标识符
    if(grepl("XENOBIOTIC.*TRANSPORTER", pathway)) {
      cleaned <- paste0(cleaned, " ", substr(pathway, 4, 6))
    }
    return(cleaned)
  }
}

# 应用映射到final_data
final_data$pathway_clean <- sapply(final_data$pathway, map_pathway_names)

# 设置通路顺序并确保唯一性
pathway_order <- selected_pathways %>%
  distinct(pathway) %>%
  pull(pathway) %>%
  sapply(map_pathway_names) %>%
  unique()

# 将pathway_clean转换为factor
final_data$pathway_clean <- factor(final_data$pathway_clean, 
                                   levels = rev(pathway_order))

# 绘制完整版本图
p <- ggplot(final_data, aes(y = factor(organ, levels = rev(unique(organ))), 
                            x = factor(pathway_clean, levels = pathway_order))) +
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
    axis.text.x = element_text(angle = 45,
                               hjust = 1,
                               size = 6),
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  labs(y = "Organ", 
       x = "Pathway",
       title = "Pathway Enrichment Z-scores Across Organs")

# 保存完整版本
ggsave("selected_pathway_zscores.pdf", p, width = 14, height = 6, dpi = 300)

# 无标签版本
p_no_labels <- ggplot(final_data, aes(y = factor(organ, levels = rev(unique(organ))), 
                                      x = factor(pathway_clean, levels = pathway_order))) +
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
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank()
  ) +
  labs(y = "", 
       x = "",
       title = "Pathway Enrichment Z-scores Across Organs")

# 保存无标签版本
ggsave("selected_pathway_zscores_no_labels.pdf", p_no_labels, width = 14, height = 4.2, dpi = 300)

# 只有x轴标签版本
p_with_x_labels <- ggplot(final_data, aes(y = factor(organ, levels = rev(unique(organ))), 
                                          x = factor(pathway_clean, levels = pathway_order))) +
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
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    axis.ticks.x = element_line(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank()
  ) +
  labs(y = "", 
       x = "",
       title = "Pathway Enrichment Z-scores Across Organs")

# 保存带X轴标签版本
ggsave("selected_pathway_zscores_with_x_labels.pdf", p_with_x_labels, width = 14, height = 5.1, dpi = 300)

# 保存数据矩阵
# Z-score矩阵
zscore_matrix <- final_data %>%
  select(organ, pathway_clean, zscore) %>%
  pivot_wider(
    names_from = pathway_clean,
    values_from = zscore
  )

# P值矩阵
pvalue_matrix <- final_data %>%
  select(organ, pathway_clean, pvalue) %>%
  pivot_wider(
    names_from = pathway_clean,
    values_from = pvalue
  )

# 保存矩阵
write.csv(zscore_matrix, "zscore_matrix.csv", row.names = FALSE)
write.csv(pvalue_matrix, "pvalue_matrix.csv", row.names = FALSE)

table(final_data$pathway_clean)
table(selected_pathways$pathway)