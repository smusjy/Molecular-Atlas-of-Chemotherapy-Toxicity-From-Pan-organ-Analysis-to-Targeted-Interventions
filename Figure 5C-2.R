library(pheatmap)
library(dplyr)
library(readxl)
library(reshape2)

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

results <- read.csv("Drug_pathway_biochemical_correlations.csv")
pathway_info <- read_excel("mice_panorgan_pathway.xlsx")
pathway_names <- read_excel("mice_panorgan_pathway_2.xlsx")

pathway_dict <- setNames(pathway_names$new_name, pathway_names$old_name)

print("检查通路映射情况：")
all_pathways <- unique(pathway_info$pathway)
mapped_pathways <- all_pathways[all_pathways %in% names(pathway_dict)]
unmapped_pathways <- all_pathways[!all_pathways %in% names(pathway_dict)]

print(paste("总通路数:", length(all_pathways)))
print(paste("能找到映射的通路数:", length(mapped_pathways)))
print(paste("需要用函数处理的通路数:", length(unmapped_pathways)))

unmapped_cleaned <- sapply(unmapped_pathways, clean_pathway_names)

final_pathway_dict <- pathway_dict
for(i in seq_along(unmapped_pathways)) {
  final_pathway_dict[unmapped_pathways[i]] <- unmapped_cleaned[i]
}

results$Drug <- gsub("Fluorouracil", "5-Fluorouracil", results$Drug)

print("Drugs in results:")
print(table(results$Drug))

print("Original unique organs:")
print(unique(pathway_info$organ))

pathway_info$organ <- gsub("Fluorouracil", "5-Fluorouracil", pathway_info$organ)
pathway_info$organ <- gsub("_vs_PB$", "", pathway_info$organ)

print("Cleaned unique organs:")
print(unique(pathway_info$organ))

drug_pathway_pairs <- split(pathway_info$pathway, pathway_info$organ)

for(drug in names(drug_pathway_pairs)) {
  cat("\n药物:", drug, "\n")
  pathways <- drug_pathway_pairs[[drug]]
  final_names <- sapply(pathways, function(x) {
    if(x %in% names(pathway_dict)) {
      pathway_dict[x]
    } else {
      clean_pathway_names(x)
    }
  })
  print(data.frame(
    Original = pathways,
    Final = final_names
  ))
}

print("Drug pathway pairs names:")
print(names(drug_pathway_pairs))

print("5-Fluorouracil in results:")
print("5-Fluorouracil" %in% unique(results$Drug))
print("5-Fluorouracil in drug_pathway_pairs:")
print("5-Fluorouracil" %in% names(drug_pathway_pairs))

filtered_results <- data.frame()

for(drug in unique(results$Drug)) {
  if(drug %in% names(drug_pathway_pairs)) {
    selected_pathways <- drug_pathway_pairs[[drug]]
    temp_results <- results %>%
      filter(Drug == drug, Pathway %in% selected_pathways)
    filtered_results <- rbind(filtered_results, temp_results)
  } else {
    print(paste("Warning: Drug", drug, "not found in pathway_pairs"))
  }
}

print("Drugs in filtered results:")
print(table(filtered_results$Drug))

if(nrow(filtered_results) > 0) {
  mat <- reshape2::dcast(filtered_results, Drug + Pathway ~ Biochemical, value.var = "Correlation")
  rownames(mat) <- paste(mat$Drug, mat$Pathway, sep="_")
  mat$Drug <- NULL
  mat$Pathway <- NULL
  
  mat_t <- t(mat)
  
  annotation_col <- data.frame(
    Drug = sub("_.*", "", colnames(mat_t)),
    row.names = colnames(mat_t)
  )
  
  drug_colors <- c(
    "#FFD33A", "#6F9AE8", "#8C9196", "#D4B17F", "#FC7444",
    "#D7E49D", "#1A7DC0", "#EE5A3E", "#467330", "#70D1F4",
    "#370332", "#0A5049", "#C50912", "#92321E", "#1A9893",
    "#FE8DBF"
  )
  names(drug_colors) <- unique(annotation_col$Drug)
  
  annotation_colors <- list(
    Drug = drug_colors
  )
  
  sig_matrix <- reshape2::dcast(filtered_results, Drug + Pathway ~ Biochemical, value.var = "P_value")
  rownames(sig_matrix) <- paste(sig_matrix$Drug, sig_matrix$Pathway, sep="_")
  sig_matrix$Drug <- NULL
  sig_matrix$Pathway <- NULL
  
  sig_marks <- matrix(
    ifelse(sig_matrix < 0.001, "***",
           ifelse(sig_matrix < 0.01, "**",
                  ifelse(sig_matrix < 0.05, "*", ""))),
    nrow = nrow(sig_matrix)
  )
  rownames(sig_marks) <- rownames(sig_matrix)
  colnames(sig_marks) <- colnames(sig_matrix)
  
  sig_marks_t <- t(sig_marks)
  
  original_pathways_col <- gsub("^.*?_", "", colnames(mat_t))
  labels_col <- sapply(original_pathways_col, function(x) {
    if(x %in% names(pathway_dict)) {
      pathway_dict[x]
    } else {
      clean_pathway_names(x)
    }
  })
  
  drug_counts <- table(annotation_col$Drug)
  gaps_col <- cumsum(drug_counts)[-length(drug_counts)]
  
  p <- pheatmap(mat_t,
                color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
                breaks = seq(-1, 1, length.out = 101),
                annotation_col = annotation_col,
                annotation_colors = annotation_colors,
                display_numbers = sig_marks_t,
                fontsize = 14,
                fontsize_row = 12,
                fontsize_col = 12,
                number_color = "black",
                angle_col = 45,
                cluster_cols = FALSE,
                cluster_rows = TRUE,
                labels_col = labels_col,
                gaps_col = gaps_col,
                main = "Pathway-Biochemical Correlation Heatmap"
  )
  
  pdf("Pathway_specific_correlation_heatmap.pdf", width = 22, height = 7)
  print(p)
  dev.off()
} else {
  print("No matching data found after filtering!")
}