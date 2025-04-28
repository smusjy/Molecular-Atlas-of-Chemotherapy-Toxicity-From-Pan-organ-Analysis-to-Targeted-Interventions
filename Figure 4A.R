library(readr)
gene_counts_results_pvalue <- read_csv("gene_counts_results_pvalue.csv")

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