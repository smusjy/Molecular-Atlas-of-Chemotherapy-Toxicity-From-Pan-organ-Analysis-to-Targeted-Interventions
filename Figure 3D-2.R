library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(ggtree)
library(aplot)
library(ape)
library(reshape2)

# Clean pathway names function
clean_pathway_names <- function(pathway) {
  lowercase_words <- c("a", "an", "the", "and", "but", "or", "for", "nor", 
                       "on", "at", "to", "by", "of", "in", "with", "within")
  
  if(!grepl("^(GO_|REACTOME_|GOBP_|GOCC_|GOMF_|HALLMARK_)", pathway)) {
    return(pathway)
  }
  
  cleaned <- gsub("^(GO_|REACTOME_|GOBP_|GOCC_|GOMF_|HALLMARK_)", "", pathway)
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

# Set working directory
setwd("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 3/Figure 3D")

# Get all csv files
files <- list.files(pattern = ".*limma.*\\.csv$")

# Create list to store all raw data
all_data_list <- list()

# Read all data and clean pathway names
for(file in files) {
  organ <- gsub(".*limma_(.+)_results\\.csv", "\\1", file)
  data <- read.csv(file, header = TRUE)
  colnames(data)[1] <- "pathway"
  data$pathway <- sapply(data$pathway, clean_pathway_names)
  all_data_list[[organ]] <- data
}

# Read immune pathway list and clean pathway names
selected_pathways <- read_excel("Immunology pathway.xlsx")
pathways_to_show <- sapply(selected_pathways$Pathway, clean_pathway_names)

# Calculate pathway z-scores and p-values
calculate_zscores_and_pvalues <- function(target_organ, all_data_list) {
  target_data <- all_data_list[[target_organ]]
  
  results <- data.frame(
    pathway = target_data$pathway,
    logFC = target_data$logFC
  )
  
  for(pathway in results$pathway) {
    all_logFCs <- numeric()
    target_logFC <- target_data$logFC[target_data$pathway == pathway]
    
    for(organ in names(all_data_list)) {
      organ_data <- all_data_list[[organ]]
      if(pathway %in% organ_data$pathway) {
        all_logFCs <- c(all_logFCs, 
                        organ_data$logFC[organ_data$pathway == pathway])
      }
    }
    
    if(length(all_logFCs) > 1) {
      mean_logFC <- mean(all_logFCs, na.rm = TRUE)
      sd_logFC <- sd(all_logFCs, na.rm = TRUE)
      
      if(!is.na(sd_logFC) && sd_logFC != 0) {
        results$zscore[results$pathway == pathway] <- 
          (target_logFC - mean_logFC) / sd_logFC
      } else {
        results$zscore[results$pathway == pathway] <- NA
      }
      
      other_logFCs <- all_logFCs[all_logFCs != target_logFC]
      
      if(length(other_logFCs) >= 2) {
        tryCatch({
          test_result <- t.test(
            other_logFCs,
            mu = target_logFC,
            alternative = "two.sided"
          )
          results$pvalue[results$pathway == pathway] <- test_result$p.value
        }, error = function(e) {
          results$pvalue[results$pathway == pathway] <- NA
        })
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

# Calculate z-scores and p-values for each organ
organ_results <- list()
final_data <- data.frame()

for(organ in names(all_data_list)) {
  results <- calculate_zscores_and_pvalues(organ, all_data_list)
  organ_results[[organ]] <- results
  
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

# Prepare clustering data
df <- final_data[,c("organ", "pathway", "zscore")]
df_matrix <- reshape2::dcast(df, organ~pathway, value.var = "zscore")
rownames(df_matrix) <- df_matrix$organ
df_matrix <- df_matrix[,-1]

# Perform hierarchical clustering
organ_clust <- hclust(dist(df_matrix))
pathway_clust <- hclust(dist(t(df_matrix)))

# 将聚类树转换为phylo对象
organ_phylo <- as.phylo(organ_clust)
pathway_phylo <- as.phylo(pathway_clust)

# Create dendrograms
p_organ <- ggtree(organ_phylo, branch.length="none") + 
  xlim(NA, 7) +
  theme_tree2()

p_pathway <- ggtree(pathway_phylo, branch.length="none") + 
  xlim(NA, 7) +
  layout_dendrogram()

# Reorder data based on clustering
final_data$organ <- factor(final_data$organ, 
                           levels = rownames(df_matrix)[organ_clust$order])
final_data$pathway <- factor(final_data$pathway, 
                             levels = colnames(df_matrix)[pathway_clust$order])

# Main plot
p_main <- ggplot(final_data, 
                 aes(x = pathway, 
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
  labs(x = "Pathway", 
       y = "Organ")

# Combine plots using aplot
final_plot <- p_main %>%
  insert_left(p_organ, width = 0.2) %>%
  insert_top(p_pathway, height = 0.2)

# Save plots
ggsave("immune_pathway_zscores_clustered.pdf", 
       final_plot, 
       width = 10, 
       height = 5.4, 
       dpi = 300)

# Version without labels
p_no_labels <- p_main +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

ggsave("immune_pathway_zscores_no_labels.pdf", 
       p_no_labels, 
       width = 12, 
       height = 3, 
       dpi = 300)

# Version with x-axis labels only
p_with_x_labels <- p_main +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    axis.title.y = element_blank()
  )

ggsave("immune_pathway_zscores_with_x_labels.pdf", 
       p_with_x_labels, 
       width = 12, 
       height = 4.3, 
       dpi = 300)

# Save clustering order information
write.csv(data.frame(
  Organ = rownames(df_matrix)[organ_clust$order],
  Cluster_Order = 1:length(organ_clust$order)
), "organ_clustering_order.csv", row.names = FALSE)

write.csv(data.frame(
  Pathway = colnames(df_matrix)[pathway_clust$order],
  Cluster_Order = 1:length(pathway_clust$order)
), "pathway_clustering_order.csv", row.names = FALSE)

# Save matrices
write.csv(df_matrix, "zscore_matrix.csv")
write.csv(reshape2::dcast(final_data, organ~pathway, value.var = "pvalue"), 
          "pvalue_matrix.csv")