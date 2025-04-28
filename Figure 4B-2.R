library(tidyverse)
library(pheatmap)

files <- list.files(pattern = ".*_limma_.*_results\\.csv$")

# 提取药物名和器官名
file_info <- str_split(files, "_") %>%
  map_df(~data.frame(drug = .[1], organ = .[3], filename = paste(., collapse="_")))

# 创建拼音到器官的映射字典
organ_dict <- c(
  "gan" = "Liver",
  "gao" = "Testis", 
  "gu" = "Bone",
  "guan" = "BloodVessel",
  "ji" = "SkeletalMuscle",
  "jie" = "Colon",
  "lung" = "Lung",
  "nao" = "Brain", 
  "pi" = "Skin",
  "qian" = "Prostate",
  "shen" = "Kidney",
  "wei" = "Stomach",
  "xin" = "Heart"
)

# 替换organ列中的拼音
file_info$organ <- organ_dict[file_info$organ]

# 查看结果
head(file_info)

# 创建组合名称
combinations <- paste(file_info$drug, file_info$organ, sep="-")

# 创建空的相关系数矩阵
cor_matrix <- matrix(NA, 
                     nrow = length(combinations), 
                     ncol = length(combinations),
                     dimnames = list(combinations, combinations))

# 两两计算组合间的相关系数
n_combinations <- length(combinations)
total_iterations <- n_combinations * n_combinations
pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)
counter <- 0

for(i in seq_along(combinations)) {
  for(j in seq_along(combinations)) {
    counter <- counter + 1
    setTxtProgressBar(pb, counter)
    
    # 获取第一个组合的信息
    idx1 <- i
    drug1 <- file_info$drug[idx1]
    organ1 <- file_info$organ[idx1]
    file1 <- files[idx1]
    
    # 获取第二个组合的信息
    idx2 <- j
    drug2 <- file_info$drug[idx2]
    organ2 <- file_info$organ[idx2]
    file2 <- files[idx2]
    
    # 读取数据
    data1 <- read.table(file1, sep=",", header=TRUE)
    data2 <- read.table(file2, sep=",", header=TRUE)
    
    data1 <- data1[,c(1:2)]
    data2 <- data2[,c(1:2)]
    
    # 重命名列名以避免混淆
    colnames(data1) <- c("X", "logFC1")
    colnames(data2) <- c("X", "logFC2")
    
    # 使用merge进行内连接，只保留交集
    merge_data <- merge(data1, data2, by = "X", all = FALSE)
    
    # 计算Spearman相关系数
    cor_matrix[i, j] <- cor(merge_data$logFC1, merge_data$logFC2, method = "spearman")
  }
}

close(pb)

write.csv(cor_matrix, "Results.csv")

# 读取数据
cor_data <- read.csv("Results.csv", row.names = 1)

# 修复列名和行名
fix_names <- function(names) {
  names <- gsub("^X", "", names)  # 移除开头的X
  names <- gsub("\\.", "-", names)  # 将点替换回横杠
  return(names)
}

# 应用修复
colnames(cor_data) <- fix_names(colnames(cor_data))
rownames(cor_data) <- fix_names(rownames(cor_data))

# 转换为矩阵
cor_matrix <- as.matrix(cor_data)

# 加载包
library(RColorBrewer)
# 获取聚类结果
hclust_result <- hclust(as.dist(1 - cor_matrix))
clusters <- cutree(hclust_result, k = 4)  # 假设切成4个聚类，可以根据实际情况调整k值

# 在每个聚类内部按器官排序
new_order <- numeric(length(clusters))
current_pos <- 1

for(i in unique(clusters)) {
  # 获取当前聚类的所有索引
  cluster_indices <- which(clusters == i)
  # 获取这些位置对应的器官
  cluster_organs <- file_info$organ[cluster_indices]
  # 按器官排序
  sorted_indices <- cluster_indices[order(cluster_organs)]
  # 填入新的顺序
  new_order[current_pos:(current_pos + length(sorted_indices) - 1)] <- sorted_indices
  current_pos <- current_pos + length(sorted_indices)
}

# 重新排序相关性矩阵
cor_matrix_ordered <- cor_matrix[new_order, new_order]

# 创建注释数据框
annotation_df <- data.frame(
  Drug = factor(file_info$drug),
  Organ = factor(file_info$organ),
  row.names = combinations
)[new_order, ]

# 手动设置颜色
drug_colors <- c(
  "#FFD33A", "#6F9AE8", "#8C9196", "#D4B17F", "#FC7444",
  "#D7E49D", "#1A7DC0", "#EE5A3E", "#467330", "#70D1F4",
  "#370332", "#0A5049", "#C50912", "#92321E", "#1A9893",
  "#FE8DBF"  # 16种颜色
)

organ_colors <- c(
  "#D4B17F", "#8C9196", "#1A9893", "#D7E49D", "#C50912",
  "#370332", "#FFD33A", "#1A7DC0", "#70D1F4", "#FC7444",
  "#467330", "#FE8DBF", "#6F9AE8"
)

# 创建注释数据框
annotation_df <- data.frame(
  Drug = factor(file_info$drug),
  Organ = factor(file_info$organ),
  row.names = combinations
)

# 创建颜色映射
annotation_colors <- list(
  Drug = setNames(drug_colors, levels(annotation_df$Drug)),
  Organ = setNames(organ_colors, levels(annotation_df$Organ))
)

# 绘制热图
pheatmap(cor_matrix_ordered,
         color = colorRampPalette(c("#204A8D", "white", "#E31A1C"))(100),
         breaks = seq(-1, 1, length.out = 101),
         cluster_rows = FALSE,  # 关闭聚类因为我们已经手动排序
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_row = annotation_df,
         annotation_col = annotation_df,
         annotation_colors = annotation_colors,
         annotation_names_row = TRUE,
         annotation_names_col = TRUE,
         annotation_legend = TRUE,
         main = "Drug-Organ Response Correlation Heatmap",
         fontsize = 8)