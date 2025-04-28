#1. 统计数量####################################################################

# 加载需要的包
library(dplyr)

# 设置工作目录
setwd("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/Protein_coding")

# 获取目录下所有文件
files <- list.files(pattern = "results.csv$")

# 创建空的结果数据框
results <- data.frame()

# 遍历每个文件
for(file in files) {
  # 读取数据
  data <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
  
  # 从文件名中提取器官和药物信息
  parts <- strsplit(file, "_")[[1]]
  drug <- parts[1]
  organ <- parts[3]
  
  # 统计显著上调的基因数量
  sig_genes <- sum(data$logFC > 2 & data$P.Value < 0.05, na.rm = TRUE)
  
  # 添加到结果数据框
  results <- rbind(results, 
                   data.frame(Organ = organ,
                              Drug = drug,
                              Sig_Genes = sig_genes))
}

# 将结果转换为宽格式
result_matrix <- results %>%
  tidyr::pivot_wider(names_from = Organ,
                     values_from = Sig_Genes) %>%
  as.data.frame()

# 查看结果
print(result_matrix)

# 保存结果
write.csv(result_matrix, "upregulated_genes_summary.csv", row.names = FALSE)


#2. 热图可视化##################################################################

library(readr)
gene_counts_results_pvalue <- read_csv("upregulated_genes_summary.csv")

head(gene_counts_results_pvalue)

library(pheatmap)

# 数据预处理并转置
matrix_data <- as.matrix(gene_counts_results_pvalue[, -1])  
rownames(matrix_data) <- gene_counts_results_pvalue$Drug    

# 创建热图
pheatmap(matrix_data,
         scale = "none",              
         cluster_rows = FALSE,        
         cluster_cols = FALSE,        
         color = colorRampPalette(c("white", "#E31A1C"))(100),  
         show_rownames = TRUE,        
         show_colnames = TRUE,        
         main = "Gene Counts Heatmap",
         fontsize = 8,                
         angle_col = 45,              
         display_numbers = TRUE,      
         number_format = "%.0f",      
         number_color = "black",      # 设置数值为纯黑色
         border_color = NA,           
         cellwidth = 20,              
         cellheight = 20,             
         width = 4,                   
         height = 4,
         legend = TRUE,
         legend_labels = "left",
         border = TRUE)